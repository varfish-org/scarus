//! Implementation of the categories for the structural variant evaluation.

pub mod common;
pub mod del;
pub mod dup;
pub mod result;

use std::{path::Path, sync::Arc};

use hgvs::data::interface::Provider;
use itertools::Itertools as _;
use prost::Message as _;

use self::common::GeneOverlap;

use super::data::clingen_dosage::Data as ClingenDosageData;
use super::data::decipher_hi::Data as DecipherHiData;
use super::data::gnomad::Data as GnomadConstraintData;
use super::data::hgnc::Data as GeneIdData;

/// Evaluator for structural variants.
pub struct Evaluator {
    /// The assembly to be used.
    #[allow(dead_code)]
    assembly: biocommons_bioutils::assemblies::Assembly,

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
}

impl Evaluator {
    /// Construct from the given paths.
    ///
    /// # Arguments
    ///
    /// * `path_tx_db` - Path to Mehari transcript database
    /// * `path_hgnc` - Path to HGNC gene identifier mapping file
    /// * `path_clingen_dosage_genes` - Path to the `ClinGen_gene_curation_list_GRCh37.tsv` file.
    /// * `path_clingen_dosage_regions` - Path to the `ClinGen_region_curation_list_GRCh37.tsv` file.
    /// * `path_decipher_hi` - Path to the `decipher_hi_prediction.tsv` file.
    /// * `path_gnomad_constraints` - Path to the `gnomad_constraints.tsv` file.
    /// * `path_clinvar_minimal` - Path to the "minimal" ClinVar RocksDB directory.
    ///
    /// # Returns
    ///
    /// A new `Evaluator`.
    ///
    /// # Errors
    ///
    /// If anything goes wrong, it returns a generic `anyhow::Error`.
    #[allow(clippy::too_many_arguments)]
    pub fn new<P1, P2, P3, P4, P5, P6, P7>(
        assembly: biocommons_bioutils::assemblies::Assembly,
        path_tx_db: P1,
        path_hgnc: P2,
        path_clingen_dosage_genes: P3,
        path_clingen_dosage_regions: P4,
        path_decipher_hi: P5,
        path_gnomad_constraints: P6,
        path_clinvar_minimal: P7,
    ) -> Result<Self, anyhow::Error>
    where
        P1: AsRef<Path>,
        P2: AsRef<Path>,
        P3: AsRef<Path>,
        P4: AsRef<Path>,
        P5: AsRef<Path>,
        P6: AsRef<Path>,
        P7: AsRef<Path>,
    {
        let provider = Self::load_provider(assembly, path_tx_db.as_ref())
            .map_err(|e| anyhow::anyhow!("failed to load transcript database: {}", e))?;
        let gene_id_data = GeneIdData::new(path_hgnc)
            .map_err(|e| anyhow::anyhow!("failed to load gene identifier data: {}", e))?;
        let clingen_dosage_data = ClingenDosageData::new(
            path_clingen_dosage_genes,
            path_clingen_dosage_regions,
            &gene_id_data,
            assembly.into(),
        )
        .map_err(|e| anyhow::anyhow!("failed to load clingen dosage data: {}", e))?;
        let decipher_hi_data = DecipherHiData::load(path_decipher_hi)
            .map_err(|e| anyhow::anyhow!("failed to load decipher hi data: {}", e))?;
        let gnomad_constraint_data =
            GnomadConstraintData::load(path_gnomad_constraints, &gene_id_data)
                .map_err(|e| anyhow::anyhow!("failed to load gnomAD constraint data: {}", e))?;
        let (clinvar_db, _) = annonars::clinvar_minimal::cli::query::open_rocksdb(
            path_clinvar_minimal,
            "clinvar",
            "meta",
        )
        .map_err(|e| anyhow::anyhow!("failed to open 'minimal' ClinVar RocksDB: {}", e))?;

        Ok(Self {
            assembly,
            chrom_to_ac: Self::chrom_to_acc(assembly, &provider),
            provider,
            gene_id_data,
            clingen_dosage_data,
            decipher_hi_data,
            gnomad_constraint_data,
            clinvar_db,
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
    ) -> Result<Vec<annonars::clinvar_minimal::pbs::Record>, anyhow::Error> {
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
                let record = annonars::clinvar_minimal::pbs::Record::decode(value)?;
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

    /// Determines whether `strucvar` overlaps with any functionally important elements.
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
    fn overlapping_elements(
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
}

#[cfg(test)]
pub mod test {
    use super::Evaluator;

    /// Fixture with the global evaluator initialize with all data paths.
    #[rstest::fixture]
    pub fn global_evaluator_37() -> Evaluator {
        Evaluator::new(
            biocommons_bioutils::assemblies::Assembly::Grch37p10,
            "tests/data/strucvars/hi_ts/txs_example_hi.bin.zst",
            "tests/data/hgnc.tsv",
            "tests/data/strucvars/ClinGen_gene_curation_list_GRCh37.tsv",
            "tests/data/strucvars/ClinGen_region_curation_list_GRCh37.tsv",
            "tests/data/strucvars/decipher_hi_prediction.tsv",
            "tests/data/strucvars/gnomad_constraints.tsv",
            "tests/data/strucvars/hi_ts/clinvar/rocksdb",
        )
        .expect("could not initialize global evaluator")
    }

    #[rstest::rstest]
    fn overlapping_elements(global_evaluator_37: super::Evaluator) {
        let res = global_evaluator_37
            .overlapping_elements("1", 8412464, 8877699)
            .expect("could not obtain overlapping elements");

        insta::assert_yaml_snapshot!(res);
    }
}
