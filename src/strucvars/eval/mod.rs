//! Implementation of the categories for the structural variant evaluation.

pub mod common;
pub mod del;
pub mod dup;
pub mod result;

use std::{path::Path, sync::Arc};

use annonars::clinvar_sv::cli::query::IntervalTrees as ClinvarSvIntervalTrees;
use annonars::common::spdi;
use annonars::functional::cli::query::IntervalTrees as FunctionalIntervalTrees;
use annonars::gnomad_pbs::exac_cnv;
use annonars::gnomad_pbs::gnomad_sv2;
use annonars::gnomad_pbs::gnomad_sv4;
use annonars::gnomad_sv::cli::query::IntervalTrees as GnomadSvIntervalTrees;
use annonars::gnomad_sv::cli::query::Record as GnomadSvRecord;
use annonars::pbs::annonars::clinvar::v1::minimal::Record as ClinvarRecord;
use annonars::pbs::annonars::clinvar::v1::sv::Record as ClinvarSvRecord;
use annonars::pbs::annonars::functional::v1::refseq::Record as FunctionalRecord;
use bio::bio_types::genome::{AbstractInterval, Interval};
use hgvs::data::interface::Provider;
use itertools::Itertools as _;
use prost::Message as _;

use self::common::GeneOverlap;

use super::data::clingen_dosage::{Data as ClingenDosageData, Gene, Region};
use super::data::decipher_hi::Data as DecipherHiData;
use super::data::gnomad::Data as GnomadConstraintData;
use super::data::hgnc::Data as GeneIdData;

/// Configuration for the evaluator.
#[derive(Debug, Clone)]
pub struct Config {
    /// Minimal reciprocal overlap for benign ClinVar SV records.
    pub clinvar_sv_min_overlap_benign: f32,
    /// Minimal reciprocal overlap for pathogenic ClinVar SV records.
    pub clinvar_sv_min_overlap_pathogenic: f32,
    /// Minimal reciprocal overlap for gnomAD SV records.
    pub gnomad_sv_min_overlap: f32,
    /// Minimal frequency in gnomAD to consider as "common".
    pub gnomad_sv_min_freq_common: f32,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            clinvar_sv_min_overlap_benign: 0.8,
            clinvar_sv_min_overlap_pathogenic: 0.8,
            gnomad_sv_min_overlap: 0.8,
            gnomad_sv_min_freq_common: 0.01,
        }
    }
}

/// Path specification for `Evaluator``.
#[derive(Debug, Clone, PartialEq, Eq, serde::Deserialize, serde::Serialize)]
pub struct Paths {
    pub path_tx_db: String,
    pub path_hgnc: String,
    pub path_clingen_dosage_genes: String,
    pub path_clingen_dosage_regions: String,
    pub path_decipher_hi: String,
    pub path_gnomad_constraints: String,
    pub path_clinvar_minimal: String,
    pub path_functional: String,
    pub path_clinvar_sv: String,
    pub path_gnomad_sv_genomes: String,
    pub path_gnomad_sv_exomes: String,
}

/// Evaluator for structural variants.
pub struct Evaluator {
    /// The assembly to be used.
    #[allow(dead_code)]
    assembly: biocommons_bioutils::assemblies::Assembly,
    /// Configuration
    config: Config,

    /// Mehari data/transcript provider.
    provider: Arc<mehari::annotate::seqvars::provider::Provider>,
    /// Mapping from chromosome name to accession for the given assembly.
    chrom_to_ac: rustc_hash::FxHashMap<String, String>,

    /// The gene identifier data.
    gene_id_data: GeneIdData,
    /// The ClinGen dosage data.
    clingen_dosage_data: ClingenDosageData,
    /// The DECIPHER HI data.
    decipher_hi_data: DecipherHiData,
    /// The gnomAD constraint data.
    gnomad_constraint_data: GnomadConstraintData,

    /// The ClinVar minimal database.
    clinvar_db: Arc<rocksdb::DBWithThreadMode<rocksdb::MultiThreaded>>,
    /// The ClinVar SV database wrapped/indexed by `IntervalTrees`.
    clinvar_sv: ClinvarSvIntervalTrees,
    /// The gnomAD SV (genomes) database wrapped/indexed by `IntervalTrees`.
    gnomad_sv_genomes: GnomadSvIntervalTrees,
    /// The gnomAD SV (exomes) database wrapped/indexed by `IntervalTrees`.
    gnomad_sv_exomes: GnomadSvIntervalTrees,
    /// The functional elements database wrapped/indexed by `IntervalTrees`.
    functional: FunctionalIntervalTrees,
}

impl Evaluator {
    /// Construct from the given paths.
    ///
    /// # Arguments
    ///
    /// * `paths` - Path specification.
    ///
    /// # Returns
    ///
    /// A new `Evaluator`.
    ///
    /// # Errors
    ///
    /// If anything goes wrong, it returns a generic `anyhow::Error`.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        assembly: biocommons_bioutils::assemblies::Assembly,
        paths: Paths,
        config: Config,
    ) -> Result<Self, anyhow::Error> {
        let provider = Self::load_provider(assembly, paths.path_tx_db.as_ref())
            .map_err(|e| anyhow::anyhow!("failed to load transcript database: {}", e))?;
        let gene_id_data = GeneIdData::new(paths.path_hgnc)
            .map_err(|e| anyhow::anyhow!("failed to load gene identifier data: {}", e))?;
        let clingen_dosage_data = ClingenDosageData::new(
            paths.path_clingen_dosage_genes,
            paths.path_clingen_dosage_regions,
            &gene_id_data,
            assembly.into(),
        )
        .map_err(|e| anyhow::anyhow!("failed to load clingen dosage data: {}", e))?;
        let decipher_hi_data = DecipherHiData::load(paths.path_decipher_hi)
            .map_err(|e| anyhow::anyhow!("failed to load decipher hi data: {}", e))?;
        let gnomad_constraint_data =
            GnomadConstraintData::load(paths.path_gnomad_constraints, &gene_id_data)
                .map_err(|e| anyhow::anyhow!("failed to load gnomAD constraint data: {}", e))?;
        let (clinvar_db, _) = annonars::clinvar_minimal::cli::query::open_rocksdb(
            &paths.path_clinvar_minimal,
            "clinvar",
            "meta",
            "clinvar_by_accession",
        )
        .map_err(|e| anyhow::anyhow!("failed to open 'minimal' ClinVar RocksDB: {}", e))?;

        let (clinvar_sv_db, clinvar_sv_meta) = annonars::clinvar_sv::cli::query::open_rocksdb(
            &paths.path_clinvar_sv,
            "clinvar_sv",
            "meta",
            "clinvar_sv_by_rcv",
        )
        .map_err(|e| anyhow::anyhow!("failed to open clinvar_sv RocksDB: {}", e))?;
        let clinvar_sv =
            ClinvarSvIntervalTrees::with_db(clinvar_sv_db, "clinvar_sv", clinvar_sv_meta)
                .map_err(|e| anyhow::anyhow!("failed to load ClinVar SV data: {}", e))?;

        let (gnomad_sv_genomes_db, gnomad_sv_genomes_meta) =
            annonars::gnomad_sv::cli::query::open_rocksdb(
                &paths.path_gnomad_sv_genomes,
                "gnomad_sv",
                "meta",
            )
            .map_err(|e| anyhow::anyhow!("failed to open gnomad_sv RocksDB: {}", e))?;
        let gnomad_sv_genomes = GnomadSvIntervalTrees::with_db(
            gnomad_sv_genomes_db,
            "gnomad_sv",
            gnomad_sv_genomes_meta,
        )
        .map_err(|e| anyhow::anyhow!("failed to load gnomAD SV data: {}", e))?;

        let (gnomad_sv_exomes_db, gnomad_sv_exomes_meta) =
            annonars::gnomad_sv::cli::query::open_rocksdb(
                &paths.path_gnomad_sv_exomes,
                "gnomad_sv",
                "meta",
            )
            .map_err(|e| anyhow::anyhow!("failed to open gnomad_sv RocksDB: {}", e))?;
        let gnomad_sv_exomes =
            GnomadSvIntervalTrees::with_db(gnomad_sv_exomes_db, "gnomad_sv", gnomad_sv_exomes_meta)
                .map_err(|e| anyhow::anyhow!("failed to load gnomAD SV data: {}", e))?;

        let (functional_db, functional_meta) = annonars::functional::cli::query::open_rocksdb(
            &paths.path_functional,
            "functional",
            "meta",
        )
        .map_err(|e| anyhow::anyhow!("failed to open functional element RocksDB: {}", e))?;
        let functional =
            FunctionalIntervalTrees::with_db(functional_db, "functional", functional_meta)
                .map_err(|e| anyhow::anyhow!("failed to load functional element data: {}", e))?;

        Ok(Self {
            config,
            assembly,
            chrom_to_ac: Self::chrom_to_acc(assembly, &provider),
            provider,
            gene_id_data,
            clingen_dosage_data,
            decipher_hi_data,
            gnomad_constraint_data,
            clinvar_db,
            clinvar_sv,
            gnomad_sv_genomes,
            gnomad_sv_exomes,
            functional,
        })
    }

    fn load_provider(
        assembly: biocommons_bioutils::assemblies::Assembly,
        path_tx_db: &Path,
    ) -> Result<Arc<mehari::annotate::seqvars::provider::Provider>, anyhow::Error> {
        tracing::info!("Opening transcript database");
        let tx_db = mehari::annotate::seqvars::load_tx_db(&format!("{}", path_tx_db.display(),))
            .map_err(|e| anyhow::anyhow!("failed to load transcript database: {}", e))?;
        tracing::info!("Building transcript interval trees ...");
        let provider = Arc::new(mehari::annotate::seqvars::provider::Provider::new(
            tx_db,
            assembly,
            mehari::annotate::seqvars::provider::ConfigBuilder::default()
                // We only consider MANE/ManePlusClinical or longest transcript.
                .transcript_picking(true)
                .build()
                .unwrap(),
        ));
        tracing::info!("... done building transcript interval trees");

        Ok(provider)
    }

    fn chrom_to_acc(
        assembly: biocommons_bioutils::assemblies::Assembly,
        provider: &Arc<mehari::annotate::seqvars::provider::Provider>,
    ) -> rustc_hash::FxHashMap<String, String> {
        let acc_to_chrom = provider.get_assembly_map(assembly);
        let mut chrom_to_acc = rustc_hash::FxHashMap::default();
        for (acc, chrom) in &acc_to_chrom {
            let chrom = if chrom.starts_with("chr") {
                chrom.strip_prefix("chr").unwrap()
            } else {
                chrom
            };
            chrom_to_acc.insert(chrom.to_string(), acc.clone());
            chrom_to_acc.insert(format!("chr{chrom}"), acc.clone());
        }

        chrom_to_acc
    }

    /// Return the gene identifier data.
    ///
    /// # Returns
    ///
    /// The gene identifier data.
    pub fn gene_id_data(&self) -> &GeneIdData {
        &self.gene_id_data
    }

    /// Query for ClinVar variants in 1-based range.
    ///
    /// # Arguments
    ///
    /// * `chrom` - Chromosome name.
    /// * `start` - Start position (1-based).
    /// * `stop` - Stop position (1-based).
    ///
    /// # Returns
    ///
    /// The ClinVar variants in the given range.
    ///
    /// # Errors
    ///
    /// If anything goes wrong, it returns a generic `anyhow::Error`.
    pub fn query_clinvar_range(
        &self,
        chrom: &str,
        start: u32,
        stop: u32,
    ) -> Result<Vec<ClinvarRecord>, anyhow::Error> {
        tracing::trace!("starting clinvar query {}:{}-{}", &chrom, start, stop);
        let start = start as i32;
        let stop = stop as i32;

        // Obtain iterator and seek to start.
        let cf_handle = self
            .clinvar_db
            .cf_handle("clinvar")
            .expect("missing clinvar column family");
        let mut iter = self.clinvar_db.raw_iterator_cf(&cf_handle);
        let pos_start = annonars::common::keys::Pos {
            chrom: chrom.to_string(),
            pos: start,
        };
        let key_start: Vec<u8> = pos_start.into();
        tracing::debug!("  seeking to key {:?}", &key_start);
        iter.seek(&key_start);

        // Cast stop to `keys::Pos`.
        let pos_stop = annonars::common::keys::Pos {
            chrom: chrom.to_string(),
            pos: stop,
        };
        let key_stop: Vec<u8> = pos_stop.into();
        tracing::debug!("  stop = {:?}", &key_stop);

        tracing::trace!("querying {}:{}-{}...", chrom, start, stop);
        let mut result = Vec::new();
        let before_query = std::time::Instant::now();
        while iter.valid() {
            if let Some(value) = iter.value() {
                // Stop if we are behind the range end already.
                let iter_key = iter.key().unwrap();
                let iter_pos: annonars::common::keys::Pos = iter_key.into();
                if iter_pos.chrom != chrom || iter_pos.pos > stop {
                    break;
                }

                // Otherwise, decode the value from the database.
                let record = ClinvarRecord::decode(value)?;
                result.push(record);

                // Proceed to the next database row.
                iter.next();
            } else {
                break;
            }
        }
        tracing::trace!(
            "... done querying {} records in {:?}",
            result.len(),
            before_query.elapsed()
        );

        Ok(result)
    }

    /// Determine overlapping genes.
    ///
    /// # Arguments
    ///
    /// * `chrom` - Chromosome name.
    /// * `start` - Start position (1-based).
    /// * `stop` - Stop position (1-based).
    ///
    /// # Returns
    ///
    /// Overlapping gene / transcript information.
    ///
    /// # Errors
    ///
    /// When the chromosome name could not be resolved or there was a problem with
    /// accessing the transcript database.
    fn overlapping_genes(
        &self,
        chrom: &str,
        start: u32,
        stop: u32,
    ) -> Result<Vec<GeneOverlap>, anyhow::Error> {
        // Map chromosome name (e.g., chr1) to chromosome accession in this assembly.
        let chrom_acc = self
            .chrom_to_ac
            .get(chrom)
            .ok_or_else(|| anyhow::anyhow!("could not resolve chromosome name `{}`", chrom))?;

        // Obtain the overlapping transcripts for the given region.
        let txs = self
            .provider
            .get_tx_for_region(chrom_acc, "splign", start as i32, stop as i32)
            .map_err(|e| anyhow::anyhow!("problem query transcript database with range: {}", e))?;

        // Extract HGNC ids / transcript ids of coding transcripts.
        let mut gene_txs = txs
            .into_iter()
            .filter_map(|tx| {
                let tx_info = self
                    .provider
                    .get_tx_info(&tx.tx_ac, &tx.alt_ac, "splign")
                    .expect("no tx info?");
                assert!(tx_info.hgnc.starts_with("HGNC:"));
                if tx_info.cds_start_i.is_some() && tx_info.cds_end_i.is_some() {
                    Some((tx_info.hgnc, tx_info.tx_ac))
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();
        gene_txs.sort();

        // Group by and collect for each gene HGNC ID.
        let gene_ovls = gene_txs
            .into_iter()
            .group_by(|(hgnc, _)| hgnc.clone())
            .into_iter()
            .map(|(hgnc, group)| {
                let txs = group.map(|(_, tx)| tx).collect::<Vec<_>>();
                let gene = self
                    .gene_id_data
                    .by_hgnc_id(&hgnc)
                    .expect("could not resolve HGNC ID")
                    .clone();
                GeneOverlap::new(gene, txs)
            })
            .collect::<Vec<_>>();

        tracing::trace!(
            "found gene overlaps for strucvar {}:{}-{}: {:?}",
            chrom,
            start,
            stop,
            &gene_ovls
        );

        Ok(gene_ovls)
    }

    /// Return (sorted) HGNC IDs of overlapping genes.
    pub fn overlapping_gene_hcnc_ids(
        &self,
        chrom: &str,
        start: u32,
        stop: u32,
    ) -> Result<Vec<String>, anyhow::Error> {
        let mut hgncs = self
            .overlapping_genes(chrom, start, stop)
            .map_err(|e| anyhow::anyhow!("failed to obtain overlapping genes: {}", e))?
            .into_iter()
            .map(|element| element.gene.hgnc_id)
            .collect::<Vec<_>>();
        hgncs.sort();
        Ok(hgncs)
    }

    /// Returns whether the given HGNC ID is protein-coding.
    ///
    /// # Arguments
    ///
    /// - `hgnc_id` -- HGNC ID to use in query.
    ///
    /// # Returns
    ///
    /// Whether the gene is protein-coding, `None`` if gene was not found.
    ///
    /// # Errors
    ///
    /// If there was an issue with querying for gene data.
    pub fn protein_coding(&self, hgnc_id: &str) -> Result<Option<bool>, anyhow::Error> {
        match self.provider.get_tx_for_gene(hgnc_id) {
            Ok(txs) => Ok(Some(txs.iter().any(|tx| tx.cds_start_i.is_some()))),
            Err(hgvs::data::error::Error::NoGeneFound(_)) => Ok(None),
            Err(e) => anyhow::bail!("problem querying transcript database: {}", e),
        }
    }

    /// Query ClinGen by overlap.
    pub fn clingen_overlaps(
        &self,
        sv_interval: &bio::bio_types::genome::Interval,
    ) -> (Vec<Gene>, Vec<Region>) {
        let mut clingen_genes = self.clingen_dosage_data.gene_by_overlap(sv_interval);
        clingen_genes.sort_by(|a, b| a.ncbi_gene_id.cmp(&b.ncbi_gene_id));
        let mut clingen_regions = self.clingen_dosage_data.region_by_overlap(sv_interval);
        clingen_regions.sort_by(|a, b| a.isca_id.cmp(&b.isca_id));
        tracing::debug!(
            "overlaps with {} ClinGen regions and {} genes",
            clingen_regions.len(),
            clingen_genes.len()
        );
        (clingen_genes, clingen_regions)
    }

    /// Query functional elements by overlap.
    pub fn functional_overlaps(
        &self,
        chrom: &str,
        start: u32,
        stop: u32,
    ) -> Result<Vec<FunctionalRecord>, anyhow::Error> {
        tracing::trace!(
            "starting functional element query {}:{}-{}",
            &chrom,
            start,
            stop
        );
        let start = start as i32;
        let stop = stop as i32;

        self.functional
            .query(&spdi::Range::new(chrom.replace("chr", ""), start, stop))
            .map_err(|e| anyhow::anyhow!("failed to query functional elements: {}", e))
    }

    /// Query ClinVar SV records by overlap.
    pub fn clinvar_sv_overlaps(
        &self,
        chrom: &str,
        start: u32,
        stop: u32,
    ) -> Result<Vec<ClinvarSvRecord>, anyhow::Error> {
        tracing::trace!("starting ClinVar SV query {}:{}-{}", &chrom, start, stop);
        let start = start as i32;
        let stop = stop as i32;

        self.clinvar_sv
            .query(&spdi::Range::new(chrom.replace("chr", ""), start, stop))
            .map_err(|e| anyhow::anyhow!("failed to query ClinVar SV db: {}", e))
    }

    /// Query gnomAD SV records by overlap.
    pub fn gnomad_sv_overlaps(
        &self,
        chrom: &str,
        start: u32,
        stop: u32,
    ) -> Result<Vec<GnomadSvRecord>, anyhow::Error> {
        tracing::trace!("starting gnomAD SV query {}:{}-{}", &chrom, start, stop);
        let start = start as i32;
        let stop = stop as i32;

        let mut res_exomes = self
            .gnomad_sv_exomes
            .query(&spdi::Range::new(chrom.replace("chr", ""), start, stop))
            .map_err(|e| anyhow::anyhow!("failed to query gnomAD SV db: {}", e))?;
        let mut res_genomes = self
            .gnomad_sv_genomes
            .query(&spdi::Range::new(chrom.replace("chr", ""), start, stop))
            .map_err(|e| anyhow::anyhow!("failed to query gnomAD SV db: {}", e))?;

        let mut result = Vec::new();
        result.append(&mut res_exomes);
        result.append(&mut res_genomes);

        Ok(result)
    }
}

pub trait IntoInterval {
    fn into_interval(self) -> Interval;
}

impl IntoInterval for &ClinvarSvRecord {
    fn into_interval(self) -> Interval {
        Interval::new(
            self.chromosome.clone(),
            (self.start as u64 - 1)..(self.stop as u64),
        )
    }
}

impl IntoInterval for &annonars::gnomad_sv::cli::query::Record {
    fn into_interval(self) -> Interval {
        match self {
            annonars::gnomad_sv::cli::query::Record::ExacCnv(record) => {
                let start = record.start as u64 - 1;
                let stop = record.stop as u64;
                let chrom = record.chrom.clone();
                Interval::new(chrom, start..stop)
            }
            annonars::gnomad_sv::cli::query::Record::GnomadSv2(record) => {
                let start = record.pos as u64 - 1;
                let stop = record.end.unwrap_or(record.pos) as u64;
                let chrom = record.chrom.clone();
                Interval::new(chrom, start..stop)
            }
            annonars::gnomad_sv::cli::query::Record::GnomadCnv4(record) => {
                let start = record.start as u64 - 1;
                let stop = record.stop as u64;
                let chrom = record.chrom.clone();
                Interval::new(chrom, start..stop)
            }
            annonars::gnomad_sv::cli::query::Record::GnomadSv4(record) => {
                let start = record.pos as u64 - 1;
                let stop = record.end.unwrap_or(record.pos) as u64;
                let chrom = record.chrom.clone();
                Interval::new(chrom, start..stop)
            }
        }
    }
}

/// Local consensus SV type.
#[derive(
    Debug,
    Clone,
    Copy,
    Default,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    serde::Deserialize,
    serde::Serialize,
)]
#[serde(rename_all = "snake_case")]
pub enum SvType {
    /// Deletion
    Del,
    /// Duplication
    Dup,
    /// Copy number variable region
    #[default]
    Cnv,
}

impl TryFrom<&annonars::gnomad_sv::cli::query::Record> for SvType {
    type Error = anyhow::Error;

    fn try_from(value: &annonars::gnomad_sv::cli::query::Record) -> Result<Self, Self::Error> {
        Ok(match value {
            GnomadSvRecord::ExacCnv(record) => {
                if record.sv_type == exac_cnv::CnvType::Del as i32 {
                    SvType::Del
                } else if record.sv_type == exac_cnv::CnvType::Dup as i32 {
                    SvType::Dup
                } else {
                    anyhow::bail!("unknown ExAC CNV type: {}", record.sv_type)
                }
            }
            GnomadSvRecord::GnomadSv2(record) => {
                if record.sv_type == gnomad_sv2::SvType::Del as i32 {
                    SvType::Del
                } else if record.sv_type == gnomad_sv2::SvType::Dup as i32 {
                    SvType::Dup
                } else if record.sv_type == gnomad_sv2::SvType::Mcnv as i32 {
                    SvType::Cnv
                } else {
                    anyhow::bail!("unable to map gnomAD v2 SV type: {}", record.sv_type)
                }
            }
            GnomadSvRecord::GnomadCnv4(record) => {
                if record.sv_type == exac_cnv::CnvType::Del as i32 {
                    SvType::Del
                } else if record.sv_type == exac_cnv::CnvType::Dup as i32 {
                    SvType::Dup
                } else {
                    anyhow::bail!("unknown gnomAD CNV type: {}", record.sv_type)
                }
            }
            GnomadSvRecord::GnomadSv4(record) => {
                if record.sv_type == gnomad_sv4::SvType::Del as i32 {
                    SvType::Del
                } else if record.sv_type == gnomad_sv4::SvType::Dup as i32 {
                    SvType::Dup
                } else if record.sv_type == gnomad_sv4::SvType::Cnv as i32 {
                    SvType::Cnv
                } else {
                    anyhow::bail!("unable to map gnomAD v4 SV type: {}", record.sv_type)
                }
            }
        })
    }
}

/// Carrier frequency in gnomAD-SV.
pub trait CarrierFrequency {
    fn carrier_freq(&self) -> Result<f32, anyhow::Error>;
}

impl CarrierFrequency for annonars::gnomad_sv::cli::query::Record {
    fn carrier_freq(&self) -> Result<f32, anyhow::Error> {
        Ok(match self {
            GnomadSvRecord::ExacCnv(_) => anyhow::bail!("ExAC CNV records must be aggregated..."),
            GnomadSvRecord::GnomadSv2(record) => {
                let allele_counts = &record
                    .allele_counts
                    .first()
                    .as_ref()
                    .expect("no allele counts")
                    .by_sex
                    .as_ref()
                    .expect("no counts by sex")
                    .overall
                    .as_ref()
                    .expect("no overall counts");
                let total = allele_counts.an - allele_counts.n_bi_genos;
                let carriers = allele_counts.n_het + allele_counts.n_homalt;
                carriers as f32 / total as f32
            }
            GnomadSvRecord::GnomadCnv4(record) => {
                let counts = record
                    .counts
                    .as_ref()
                    .expect("no counts")
                    .overall
                    .as_ref()
                    .expect("no counts by sex");
                counts.sc as f32 / counts.sn as f32
            }
            GnomadSvRecord::GnomadSv4(record) => {
                let allele_counts = record
                    .allele_counts
                    .first()
                    .as_ref()
                    .expect("no allele counts")
                    .by_sex
                    .as_ref()
                    .expect("no counts by sex")
                    .overall
                    .as_ref()
                    .expect("no overall counts");
                let total = allele_counts.an - allele_counts.n_bi_genos;
                let carriers = allele_counts.n_het + allele_counts.n_homalt;
                carriers as f32 / total as f32
            }
        })
    }
}

/// Reciprocal overlap between two intervals.
pub trait Overlaps {
    fn overlap_len(&self, other: &Self) -> u64;
    fn reciprocal_overlap(&self, other: &Self) -> f64;
}

impl Overlaps for Interval {
    fn overlap_len(&self, other: &Self) -> u64 {
        if self.contig() != other.contig() {
            0
        } else {
            let start = std::cmp::max(self.range().start, other.range().start);
            let end = std::cmp::min(self.range().end, other.range().end);
            if start < end {
                end - start
            } else {
                0
            }
        }
    }

    fn reciprocal_overlap(&self, other: &Self) -> f64 {
        let self_len = self.range().end - self.range().start;
        let other_len = other.range().end - other.range().start;
        let overlap_len = self.overlap_len(other) as f64;
        overlap_len / std::cmp::max(self_len, other_len) as f64
    }
}

#[cfg(test)]
pub mod test {
    use crate::strucvars::eval::Overlaps;

    use super::Evaluator;

    #[rstest::rstest]
    #[case("1", 1, 10, "1", 5, 15, 5)]
    #[case("1", 1, 10, "1", 10, 15, 0)]
    #[case("1", 1, 10, "2", 1, 10, 0)]
    fn overlap_len(
        #[case] contig1: &str,
        #[case] start1: u64,
        #[case] stop1: u64,
        #[case] contig2: &str,
        #[case] start2: u64,
        #[case] stop2: u64,
        #[case] expected: u64,
    ) {
        assert_eq!(
            super::Interval::new(contig1.into(), start1..stop1)
                .overlap_len(&super::Interval::new(contig2.into(), start2..stop2)),
            expected
        );
    }

    #[rstest::rstest]
    #[case("1", 0, 10, "1", 5, 15, 0.5)]
    #[case("1", 0, 10, "1", 10, 20, 0.0)]
    #[case("1", 0, 10, "2", 0, 10, 0.0)]
    #[case("1", 0, 10, "1", 5, 25, 0.25)]
    fn reciprocal_overlap(
        #[case] contig1: &str,
        #[case] start1: u64,
        #[case] stop1: u64,
        #[case] contig2: &str,
        #[case] start2: u64,
        #[case] stop2: u64,
        #[case] expected: f64,
    ) {
        assert_eq!(
            super::Interval::new(contig1.into(), start1..stop1)
                .reciprocal_overlap(&super::Interval::new(contig2.into(), start2..stop2)),
            expected
        );
    }
    /// Fixture with the global evaluator initialize with all data paths.
    #[rstest::fixture]
    #[once]
    pub fn global_evaluator_37() -> Evaluator {
        Evaluator::new(
            biocommons_bioutils::assemblies::Assembly::Grch37p10,
            super::Paths {
                path_tx_db: "tests/data/strucvars/hi_ts/txs_example_hi.bin.zst".into(),
                path_hgnc: "tests/data/hgnc.tsv".into(),
                path_clingen_dosage_genes:
                    "tests/data/strucvars/ClinGen_gene_curation_list_GRCh37.tsv".into(),
                path_clingen_dosage_regions:
                    "tests/data/strucvars/ClinGen_region_curation_list_GRCh37.tsv".into(),
                path_decipher_hi: "tests/data/strucvars/decipher_hi_prediction.tsv".into(),
                path_gnomad_constraints: "tests/data/strucvars/gnomad_constraints.tsv".into(),
                path_clinvar_minimal: "tests/data/strucvars/hi_ts/clinvar/rocksdb".into(),
                path_functional: "tests/data/strucvars/hi_ts/functional/rocksdb".into(),
                path_clinvar_sv: "tests/data/strucvars/hi_ts/clinvar-sv/rocksdb".into(),
                path_gnomad_sv_genomes: "tests/data/strucvars/hi_ts/gnomad-sv/gnomad-sv2/rocksdb".into(),
                path_gnomad_sv_exomes: "tests/data/strucvars/hi_ts/gnomad-sv/exac-cnv/rocksdb".into(),
            },
            Default::default(),
        )
        .expect("could not initialize global evaluator")
    }

    #[rstest::rstest]
    fn overlapping_elements(global_evaluator_37: &super::Evaluator) {
        let res = global_evaluator_37
            .overlapping_genes("1", 8412464, 8877699)
            .expect("could not obtain overlapping elements");

        insta::assert_yaml_snapshot!(res);
    }

    /// Test internal working of `clingen_overlaps`.
    #[tracing_test::traced_test]
    #[rstest::rstest]
    #[case("17", 67_892_000, 69_793_000, "region-match")]
    #[case("X", 152_990_000, 153_011_000, "gene-abcd1")]
    fn clingen_overlaps(
        #[case] chrom: &str,
        #[case] start: u64,
        #[case] stop: u64,
        #[case] label: &str,
        global_evaluator_37: &Evaluator,
    ) -> Result<(), anyhow::Error> {
        mehari::common::set_snapshot_suffix!("{}", label);

        let sv_interval =
            crate::strucvars::data::intervals::Interval::new(chrom.into(), start..stop);

        let (genes, regions) = global_evaluator_37.clingen_overlaps(&sv_interval);

        insta::assert_yaml_snapshot!(genes);
        insta::assert_yaml_snapshot!(regions);

        Ok(())
    }

    /// Test internal working of `functional_overlaps`.
    #[tracing_test::traced_test]
    #[rstest::rstest]
    #[case("1", 1_004_592, 1_004_706, "chr1-example")]
    #[case("20", 3_682_289, 6_749_387, "chr20-example")]
    fn functional_overlaps(
        #[case] chrom: &str,
        #[case] start: u32,
        #[case] stop: u32,
        #[case] label: &str,
        global_evaluator_37: &Evaluator,
    ) -> Result<(), anyhow::Error> {
        mehari::common::set_snapshot_suffix!("{}", label);

        let result = global_evaluator_37.functional_overlaps(chrom, start, stop)?;
        insta::assert_yaml_snapshot!(result);

        Ok(())
    }

    /// Test internal working of `clinvar_sv_overlaps`.
    #[tracing_test::traced_test]
    #[rstest::rstest]
    #[case("1", 8_018_845, 8_044_954, "VCV000007063")]
    #[case("21", 1, 1, "almost-empty")]
    fn clinvar_sv_overlaps(
        #[case] chrom: &str,
        #[case] start: u32,
        #[case] stop: u32,
        #[case] label: &str,
        global_evaluator_37: &Evaluator,
    ) -> Result<(), anyhow::Error> {
        mehari::common::set_snapshot_suffix!("{}", label);

        let result = global_evaluator_37.clinvar_sv_overlaps(chrom, start, stop)?;
        insta::assert_yaml_snapshot!(result);

        Ok(())
    }

    /// Test internal working of `gnomad_sv_overlap`.
    #[tracing_test::traced_test]
    #[rstest::rstest]
    // DELs
    //
    // ExAC-CNV
    #[case("10", 51_620_321, 51_748_701, "chr10-51620321-51748701-NFE")]
    // gnomAD-SV v2
    #[case("1", 21000, 26000, "gnomAD-SV_v2.1_DEL_1_")]
    // DUPs
    //
    // ExAC-CNV
    #[case("10", 92_994, 95_183, "chr10-92994-95183-NFE")]
    // gnomAD-SV v2
    #[case("1", 157_000, 166_000, "gnomAD-SV_v2.1_DUP_1_6")]
    fn gnomad_sv_overlaps(
        #[case] chrom: &str,
        #[case] start: u32,
        #[case] stop: u32,
        #[case] label: &str,
        global_evaluator_37: &Evaluator,
    ) -> Result<(), anyhow::Error> {
        mehari::common::set_snapshot_suffix!("{}", label);

        let result = global_evaluator_37.gnomad_sv_overlaps(chrom, start, stop)?;
        insta::assert_yaml_snapshot!(result);

        Ok(())
    }

    /// Test `is_protein_coding`.
    #[tracing_test::traced_test]
    #[rstest::rstest]
    // currently RNA genes not imported
    // #[case("HGNC:53963", Some(false))]
    #[case("HGNC:39941", Some(true))]
    #[case("HGNC:xxx", None)]
    fn is_protein_coding(
        #[case] hgnc_id: &str,
        #[case] expected: Option<bool>,
        global_evaluator_37: &super::Evaluator,
    ) -> Result<(), anyhow::Error> {
        mehari::common::set_snapshot_suffix!("{}", hgnc_id);

        assert_eq!(global_evaluator_37.protein_coding(hgnc_id)?, expected);

        Ok(())
    }
}
