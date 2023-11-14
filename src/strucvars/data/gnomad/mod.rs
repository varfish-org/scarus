//! gnomAD gene constraints.

//! Code for accessing gnomAD gene constraints data.

use std::path::Path;

use super::hgnc;

pub mod io;

/// The HI informative part of gnomAD constraints.
#[serde_with::skip_serializing_none]
#[derive(Debug, Clone, serde::Deserialize, serde::Serialize)]
pub struct Record {
    /// The HGNC gene identifier.
    pub hgnc_id: String,
    /// pLI score.
    pub pli: Option<f64>,
    /// The upper bound of the loss-of-function observed/expected CI.
    pub oe_lof_upper: Option<f64>,
}

impl Record {
    /// Create a new `Record` from `io::Record`.
    ///
    /// # Arguments
    ///
    /// * `other` - The `io::Record` to convert.
    /// * `hgnc_id` - The HGNC gene identifier.
    ///
    /// # Returns
    ///
    /// A new `Record`.
    pub fn from_io_record(other: io::Record, hgnc_id: String) -> Self {
        Self {
            hgnc_id,
            pli: other.pli,
            oe_lof_upper: other.oe_lof_upper,
        }
    }
}

/// Load gnomAD constraints file and map ENSEMBL gene identifiers.
///
/// # Arguments
///
/// * `path` - Path to gnomAD constraints file.
/// * `hgnc_data` - `hgnc::Data` to use for gene ID mapping.
///
/// # Returns
///
/// gnomAD constraints.
///
/// # Errors
///
/// If anything goes wrong, it returns a generic `anyhow::Error`.
pub fn load_file<P>(path: P, hgnc_data: &hgnc::Data) -> Result<Vec<Record>, anyhow::Error>
where
    P: AsRef<Path>,
{
    Ok(io::load_file(path)?
        .into_iter()
        .flat_map(|record| {
            hgnc_data
                .by_ensembl_gene_id(&record.ensembl_gene_id)
                .map(|gene_id_info| Record::from_io_record(record, gene_id_info.hgnc_id.clone()))
        })
        .collect::<Vec<_>>())
}

/// Facade struct that allows easy acccess to the gnomAD constraint data.
#[derive(Debug, Clone)]
pub struct Data {
    /// The gnomAD constraint.
    data: Vec<Record>,
    /// Mapping from HGNC identifier to gnomAD constraint data index.
    hgnc_to_data_idx: rustc_hash::FxHashMap<String, usize>,
}

impl Data {
    /// Load from file and construct.
    pub fn load<P>(path: P, hgnc_data: &hgnc::Data) -> Result<Self, anyhow::Error>
    where
        P: AsRef<Path>,
    {
        let data = load_file(path, hgnc_data)?;
        Ok(Self::new(data))
    }

    /// Create a new `Data` object.
    ///
    /// # Arguments
    ///
    /// * `data` - The gnomAD constraint record.
    ///
    /// # Returns
    ///
    /// A new `Data` object.
    pub fn new(data: Vec<Record>) -> Self {
        let hgnc_to_data_idx = data
            .iter()
            .enumerate()
            .map(|(idx, record)| (record.hgnc_id.clone(), idx))
            .collect::<rustc_hash::FxHashMap<_, _>>();
        Self {
            data,
            hgnc_to_data_idx,
        }
    }

    /// Get the gnomAD constraint record for the given HGNC identifier.
    ///
    /// # Arguments
    ///
    /// * `hgnc_id` - HGNC identifier.
    ///
    /// # Returns
    ///
    /// The gnomAD constraint record for the given HGNC identifier, if any.
    pub fn by_hgnc_id(&self, hgnc_id: &str) -> Option<&Record> {
        self.hgnc_to_data_idx
            .get(hgnc_id)
            .map(|idx| &self.data[*idx])
    }
}

#[cfg(test)]
mod test {
    #[tracing_test::traced_test]
    #[test]
    fn test_load_file() -> Result<(), anyhow::Error> {
        let records = super::load_file(
            "tests/data/strucvars/gnomad_constraints.tsv",
            &super::hgnc::Data::new("tests/data/hgnc.tsv")?,
        )?;

        assert_eq!(records.len(), 18134);
        insta::assert_yaml_snapshot!(&records[0..5]);

        Ok(())
    }

    #[tracing_test::traced_test]
    #[test]
    fn data_load() -> Result<(), anyhow::Error> {
        let data = super::Data::load(
            "tests/data/strucvars/gnomad_constraints.tsv",
            &super::hgnc::Data::new("tests/data/hgnc.tsv")?,
        )?;

        insta::assert_yaml_snapshot!(data.by_hgnc_id("HGNC:5"));

        Ok(())
    }
}
