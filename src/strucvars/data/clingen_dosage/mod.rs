//! Data structures and I/O code for `ClinGen` dosage sensitivity.

use std::path::Path;

use bio::{
    bio_types::genome::AbstractInterval as _,
    data_structures::interval_tree::ArrayBackedIntervalTree,
};
use itertools::Itertools as _;

use crate::common::Assembly;

use super::{
    hgnc,
    intervals::{do_overlap, Interval},
};

pub mod io;

/// Haploinsufficiency scores.
#[derive(
    Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, serde::Deserialize, serde::Serialize,
)]
#[serde(rename_all = "snake_case")]
pub enum Score {
    /// Sufficient evidence for dosage pathogenicity (3)
    SufficientEvidence,
    /// Little evidence for dosage pathogenicity (1)
    LittleEvidence,
    /// Some evidence for dosage pathogenicity (2)
    SomeEvidence,
    /// No evidence for dosage pathogenicity (0)
    NoEvidenceAvailable,
    /// Gene associated with autosomal recessive phenotype (30)
    GeneAssociatedWithRecessivePhenotype,
    /// Dosage sensitivity unlikely (40)
    DosageSensitivityUnlikely,
}

impl TryFrom<Option<u32>> for Score {
    type Error = anyhow::Error;

    fn try_from(value: Option<u32>) -> Result<Self, Self::Error> {
        match value {
            None | Some(0) => Ok(Self::NoEvidenceAvailable),
            Some(1) => Ok(Self::LittleEvidence),
            Some(2) => Ok(Self::SomeEvidence),
            Some(3) => Ok(Self::SufficientEvidence),
            Some(30) => Ok(Self::GeneAssociatedWithRecessivePhenotype),
            Some(40) => Ok(Self::DosageSensitivityUnlikely),
            _ => anyhow::bail!("invalid score: {:?}", value),
        }
    }
}

/// `ClinGen` dosage sensitivy gene record to be used in the app.
#[derive(Debug, Clone, PartialEq, serde::Deserialize, serde::Serialize)]
pub struct Gene {
    /// Gene symbol.
    pub gene_symbol: String,
    /// NCBI gene ID.
    pub ncbi_gene_id: String,
    /// Genomic location.
    pub genomic_location: String,
    /// Haploinsufficiency score.
    pub haploinsufficiency_score: Score,
    /// Triplosensitivity score.
    pub triplosensitivity_score: Score,
    /// Haploinsufficiency Disease ID.
    pub haploinsufficiency_disease_id: Option<String>,
    /// Haploinsufficiency Disease ID.
    pub triplosensitivity_disease_id: Option<String>,
}

/// Helper to convert genomic location string into an interval.
fn genomic_location_to_interval(
    genomic_location: &str,
) -> Result<bio::bio_types::genome::Interval, anyhow::Error> {
    let mut parts = genomic_location.split(':');
    let chrom = parts.next().ok_or_else(|| {
        anyhow::anyhow!(
            "could not parse chromosome from genomic location: {}",
            genomic_location
        )
    })?;
    let mut parts = parts
        .next()
        .ok_or_else(|| anyhow::anyhow!("could not parse region {}", genomic_location))?
        .split('-');
    let begin = parts
        .next()
        .unwrap()
        .parse::<u64>()
        .map_err(|e| anyhow::anyhow!("could not parse start position from: {}", e))?
        .saturating_sub(1);
    let end = parts
        .next()
        .unwrap()
        .parse::<u64>()
        .map_err(|e| anyhow::anyhow!("could not parse end position from: {}", e))?;
    Ok(bio::bio_types::genome::Interval::new(
        chrom.to_string(),
        begin..end,
    ))
}

impl TryInto<bio::bio_types::genome::Interval> for Gene {
    type Error = anyhow::Error;

    fn try_into(self) -> Result<bio::bio_types::genome::Interval, Self::Error> {
        genomic_location_to_interval(&self.genomic_location)
    }
}

/// `ClinGen` dosage sensitivy region record to be used in the app.
#[derive(Debug, Clone, PartialEq, serde::Deserialize, serde::Serialize)]
pub struct Region {
    /// ISCA ID
    pub isca_id: String,
    /// ISCA Region Name
    pub isca_region_name: String,
    /// Genomic locaion.
    pub genomic_location: String,
    /// Haploinsufficiency score.
    pub haploinsufficiency_score: Score,
    /// Triplosensitivity score.
    pub triplosensitivity_score: Score,
    /// Haploinsufficiency Disease ID.
    pub haploinsufficiency_disease_id: Option<String>,
    /// Haploinsufficiency Disease ID.
    pub triplosensitivity_disease_id: Option<String>,
}

impl TryInto<bio::bio_types::genome::Interval> for Region {
    type Error = anyhow::Error;

    fn try_into(self) -> Result<bio::bio_types::genome::Interval, Self::Error> {
        genomic_location_to_interval(&self.genomic_location)
    }
}

/// Facade struct for querying the `ClinGen` dosage information.
#[derive(Debug, Clone)]
pub struct Data {
    /// Information regarding genes.
    genes: Vec<Gene>,
    /// Mapping from HGNC gene ID to index in `genes`.
    hgnc_to_genes_idx: rustc_hash::FxHashMap<String, usize>,
    /// Interval tree index for querying genes, one per chromosome.
    genes_trees: rustc_hash::FxHashMap<String, ArrayBackedIntervalTree<u64, usize>>,

    /// Information on regions.
    regions: Vec<Region>,
}

impl Data {
    /// Construct `Info` from paths to `ClinGen` dosage files.
    ///
    /// # Arguments
    ///
    /// * `gene_path` - Path to the `ClinGen_gene_curation_list_GRCh37.tsv` file.
    /// * `region_path` - Path to the `ClinGen_region_curation_list_GRCh37.tsv` file.
    /// * `hgnc_data` - `hgnc::Data` to build mapping from HGNC ID to NCBI gene ID.
    /// * `assembly` - Assembly that will be used (for adjusting chromosome names).
    ///
    /// # Returns
    ///
    /// `Info` struct.
    ///
    /// # Errors
    ///
    /// If anything goes wrong, it returns a generic `anyhow::Error`.
    pub fn new<P1, P2>(
        gene_path: P1,
        region_path: P2,
        hgnc_data: &hgnc::Data,
        assembly: Assembly,
    ) -> Result<Self, anyhow::Error>
    where
        P1: AsRef<Path>,
        P2: AsRef<Path>,
    {
        let genes: Vec<Gene> = io::load_file::<io::Gene, _>(gene_path, assembly)
            .map_err(|e| anyhow::anyhow!("problem loading genes: {}", e))?
            .into_iter()
            .filter(|gene| match gene.genomic_location.as_str() {
                "tbd" => {
                    tracing::warn!(
                        "skipping gene with genomic location {}: {}",
                        gene.genomic_location,
                        gene.gene_symbol
                    );
                    false
                }
                _ => true,
            })
            .map(std::convert::TryInto::try_into)
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e| anyhow::anyhow!("problem converting genes: {}", e))?;
        let genes_trees = Self::build_genes_trees(&genes)
            .map_err(|e| anyhow::anyhow!("problem building interval trees for genes: {}", e))?;
        let hgnc_to_genes_idx = genes
            .iter()
            .enumerate()
            .map(|(idx, gene)| -> Result<(String, usize), anyhow::Error> {
                Ok((
                    hgnc_data
                        .by_ncbi_gene_id(&gene.ncbi_gene_id)
                        .map(|info| info.hgnc_id.clone())
                        .ok_or_else(|| {
                            anyhow::anyhow!(
                                "could not find HGNC ID for NCBI gene ID {}",
                                gene.ncbi_gene_id
                            )
                        })?,
                    idx,
                ))
            })
            .collect::<Result<rustc_hash::FxHashMap<_, _>, _>>()?;

        let regions = io::load_file::<io::Region, _>(region_path, assembly)
            .map_err(|e| anyhow::anyhow!("problem loading regions: {}", e))?
            .into_iter()
            .filter(|gene| match gene.genomic_location.as_str() {
                "tbd" => {
                    tracing::warn!(
                        "skipping region with genomic location {}: {}",
                        gene.genomic_location,
                        gene.isca_id
                    );
                    false
                }
                _ => true,
            })
            .map(std::convert::TryInto::try_into)
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e| anyhow::anyhow!("problem converting regions: {}", e))?;

        Ok(Self {
            genes,
            genes_trees,
            hgnc_to_genes_idx,
            regions,
        })
    }

    /// Build interval trees for the genes.
    pub fn build_genes_trees(
        genes: &[Gene],
    ) -> Result<rustc_hash::FxHashMap<String, ArrayBackedIntervalTree<u64, usize>>, anyhow::Error>
    {
        let mut per_contig = genes
            .iter()
            .enumerate()
            .map(|(idx, gene)| {
                let interval: Interval = gene.clone().try_into()?;
                let contig = interval.contig().to_string();
                let range = interval.range();
                Ok((contig, (range, idx)))
            })
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e: anyhow::Error| {
                anyhow::anyhow!("problem converting genes to intervals: {}", e)
            })?;
        per_contig.sort_by_key(|(contig, _)| contig.clone());
        let per_contig = per_contig
            .into_iter()
            .group_by(|(contig, _)| contig.clone())
            .into_iter()
            .map(|(contig, group)| {
                (
                    contig,
                    ArrayBackedIntervalTree::from_iter(group.map(|(_, (range, idx))| (range, idx))),
                )
            })
            .collect::<Vec<_>>();

        Ok(rustc_hash::FxHashMap::from_iter(per_contig))
    }

    /// Obtain the gene information for the given gene HGNC ID.
    ///
    /// # Arguments
    ///
    /// * `hgnc_id` - HGNC identifier.
    ///
    /// # Returns
    ///
    /// The ClinGen gene dosage information for the given `hgnc_id`, if any.
    pub fn gene_by_hgnc_id(&self, hgnc_id: &str) -> Option<&Gene> {
        self.hgnc_to_genes_idx
            .get(hgnc_id)
            .map(|idx| &self.genes[*idx])
    }

    /// Obtain overlapping gene information for the given genomic location.
    ///
    /// # Arguments
    ///
    /// * `locus` - Genomic location.
    ///
    /// # Returns
    ///
    /// The ClinGen gene dosage information overlapping with `locus`.
    pub fn gene_by_overlap(&self, locus: &bio::bio_types::genome::Interval) -> Vec<Gene> {
        self.genes_trees
            .get(locus.contig())
            .map(|tree| {
                tree.find(locus.range())
                    .into_iter()
                    .map(|idx| &self.genes[*idx.data()])
                    .cloned()
                    .collect::<Vec<_>>()
            })
            .unwrap_or_default()
    }

    /// Obtain overlapping region information for the given genomic location.
    ///
    /// # Arguments
    ///
    /// * `locus` - Genomic location.
    ///
    /// # Returns
    ///
    /// The ClinGen region dosage information overlapping with `locus`.
    pub fn region_by_overlap(&self, locus: &bio::bio_types::genome::Interval) -> Vec<Region> {
        self.regions
            .iter()
            .filter(|region| {
                let region_locus: bio::bio_types::genome::Interval = (*region)
                    .clone()
                    .try_into()
                    .expect("region to interval conversion must succeed");
                do_overlap(locus, &region_locus)
            })
            .cloned()
            .collect()
    }
}

#[cfg(test)]
mod test {
    use crate::common::Assembly;

    #[test]
    fn info_with_paths() -> Result<(), anyhow::Error> {
        let info = super::Data::new(
            "tests/data/strucvars/ClinGen_gene_curation_list_GRCh37.tsv",
            "tests/data/strucvars/ClinGen_region_curation_list_GRCh37.tsv",
            &super::hgnc::Data::new("tests/data/hgnc.tsv")?,
            Assembly::Grch37,
        )?;

        assert_eq!(info.genes.len(), 1518);
        assert_eq!(info.regions.len(), 513);

        Ok(())
    }

    #[test]
    fn score_conversion() -> Result<(), anyhow::Error> {
        assert_eq!(
            super::Score::try_from(Some(0))?,
            super::Score::NoEvidenceAvailable
        );
        assert_eq!(
            super::Score::try_from(Some(1))?,
            super::Score::LittleEvidence
        );
        assert_eq!(super::Score::try_from(Some(2))?, super::Score::SomeEvidence);
        assert_eq!(
            super::Score::try_from(Some(3))?,
            super::Score::SufficientEvidence
        );
        assert_eq!(
            super::Score::try_from(Some(30))?,
            super::Score::GeneAssociatedWithRecessivePhenotype
        );
        assert_eq!(
            super::Score::try_from(Some(40))?,
            super::Score::DosageSensitivityUnlikely
        );
        assert_eq!(
            super::Score::try_from(None)?,
            super::Score::NoEvidenceAvailable
        );
        assert!(super::Score::try_from(Some(4)).is_err());

        Ok(())
    }
}
