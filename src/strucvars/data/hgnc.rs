//! HGNC-based gene identifier mapping.

use std::{io::BufReader, path::Path};

/// Mapping of gene identifiers.
#[derive(Debug, Clone, Default, PartialEq, serde::Deserialize, serde::Serialize)]
pub struct GeneIdInfo {
    /// HGNC identifier.
    pub hgnc_id: String,
    /// Official HGNC gene symbol.
    #[serde(alias = "gene_symbol")]
    pub hgnc_symbol: Option<String>,
    /// ENSEMBL gene identifier.
    pub ensembl_gene_id: Option<String>,
    /// NCBI gene identifier.
    #[serde(alias = "entrez_id")]
    pub ncbi_gene_id: Option<String>,
}

/// Load `HGVS` gene mapping file.
///
/// # Arguments
///
/// * `path` - Path to the mapping TSV file with the fields
///   `hgnc_id`, `ensembl_gene_id`, `entrez_id`, `gene_symbol`.
///
/// # Returns
///
/// Mapping records.
///
/// # Errors
///
/// If anything goes wrong, it returns a generic `anyhow::Error`.
pub fn load_file<P>(path: P) -> Result<Vec<GeneIdInfo>, anyhow::Error>
where
    P: AsRef<Path>,
{
    // Construct buffered file and CSV reader.
    let reader = std::fs::File::open(path)
        .map_err(|e| anyhow::anyhow!("problem opening file: {}", e))
        .map(BufReader::new)?;
    let mut csv_reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .flexible(false)
        .from_reader(reader);
    let mut result = Vec::new();
    for record in csv_reader.deserialize() {
        let record = record.map_err(|e| anyhow::anyhow!("problem parsing record: {}", e))?;
        result.push(record);
    }

    Ok(result)
}

/// Facade struct for accessing gene identifier information.
pub struct Data {
    /// Gene identifier information.
    gene_id_infos: Vec<GeneIdInfo>,
    /// Mapping from HGNC identifier to gene identifier information.
    hgnc_to_infos_idx: rustc_hash::FxHashMap<String, usize>,
    /// Mapping from NCBI gene ID to gene identifier information.
    ncbi_to_infos_idx: rustc_hash::FxHashMap<String, usize>,
    /// Mapping from ENSEMBL gene ID to gene identifier information.
    ensg_to_infos_idx: rustc_hash::FxHashMap<String, usize>,
}

impl Data {
    /// Construct from the given path.
    ///
    /// # Arguments
    ///
    /// * `path_hgnc` - Path to HGNC gene identifier mapping file.
    ///
    /// # Returns
    ///
    /// A new `GeneIdData`.
    ///
    /// # Errors
    ///
    /// If anything goes wrong, it returns a generic `anyhow::Error`.
    pub fn new<P>(path_hgnc: P) -> Result<Self, anyhow::Error>
    where
        P: AsRef<Path>,
    {
        let gene_id_infos = load_file(path_hgnc.as_ref())?;
        let hgnc_to_infos_idx = gene_id_infos
            .iter()
            .enumerate()
            .map(|(idx, info)| (info.hgnc_id.clone(), idx))
            .collect();
        let ncbi_to_infos_idx = gene_id_infos
            .iter()
            .enumerate()
            .filter_map(|(idx, info)| info.ncbi_gene_id.as_ref().map(|id| (id.clone(), idx)))
            .collect();
        let ensg_to_infos_idx = gene_id_infos
            .iter()
            .enumerate()
            .filter_map(|(idx, info)| info.ensembl_gene_id.as_ref().map(|id| (id.clone(), idx)))
            .collect();
        Ok(Self {
            gene_id_infos,
            hgnc_to_infos_idx,
            ncbi_to_infos_idx,
            ensg_to_infos_idx,
        })
    }

    /// Obtain the gene identifier information for the given HGNC identifier.
    ///
    /// # Arguments
    ///
    /// * `hgnc_id` - HGNC identifier.
    ///
    /// # Returns
    ///
    /// The gene identifier information, if any available for `hgnc_id`.
    pub fn by_hgnc_id(&self, hgnc_id: &str) -> Option<&GeneIdInfo> {
        self.hgnc_to_infos_idx
            .get(hgnc_id)
            .map(|idx| &self.gene_id_infos[*idx])
    }

    /// Obtain the gene identifier information for the given NCBI gene ID.
    ///
    /// # Arguments
    ///
    /// * `ncbi_gene_id` - NCBI gene identifier.
    ///
    /// # Returns
    ///
    /// The gene identifier information, if any available for `hgnc_id`.
    pub fn by_ncbi_gene_id(&self, ncbi_gene_id: &str) -> Option<&GeneIdInfo> {
        self.ncbi_to_infos_idx
            .get(ncbi_gene_id)
            .map(|idx| &self.gene_id_infos[*idx])
    }

    /// Obtain the gene identifier information for the given NCBI gene ID.
    ///
    /// # Arguments
    ///
    /// * `ensembl_gene_id` - ENSEMBL gene identifier.
    ///
    /// # Returns
    ///
    /// The gene identifier information, if any available for `hgnc_id`.
    pub fn by_ensembl_gene_id(&self, ensembl_gene_id: &str) -> Option<&GeneIdInfo> {
        self.ensg_to_infos_idx
            .get(ensembl_gene_id)
            .map(|idx| &self.gene_id_infos[*idx])
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_load_file() -> Result<(), anyhow::Error> {
        let path = "tests/data/hgnc.tsv";
        let records = load_file(path)?;

        assert_eq!(records.len(), 43626);
        insta::assert_yaml_snapshot!(&records[0..5]);

        Ok(())
    }

    #[test]
    fn test_gene_id_data() -> Result<(), anyhow::Error> {
        let data = Data::new("tests/data/hgnc.tsv")?;

        assert_eq!(data.gene_id_infos.len(), 43626);
        insta::assert_yaml_snapshot!(data.by_hgnc_id("HGNC:1100"));

        Ok(())
    }
}
