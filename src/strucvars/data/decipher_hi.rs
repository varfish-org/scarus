//! Haploinsufficiency index from DECIPHER.

use std::{io::BufReader, path::Path};

/// DECIPHER HI prediction.
#[derive(Debug, Clone, PartialEq, serde::Deserialize, serde::Serialize)]
pub struct HiPrediction {
    /// HGNC identifier.
    pub hgnc_id: String,
    /// Official HGNC gene symbol.
    pub hgnc_symbol: String,
    /// P(HI) prediction from DECIPHER HI.
    pub p_hi: f64,
    /// Percent HI index.
    pub hi_index: f64,
}

/// Load DECIPHER HI regions file.
///
/// # Arguments
///
/// * `path` - Path to DECIPHER HI file.
///
/// # Returns
///
/// HI predicitons.
///
/// # Errors
///
/// If anything goes wrong, it returns a generic `anyhow::Error`.
pub fn load_file<P>(path: P) -> Result<Vec<HiPrediction>, anyhow::Error>
where
    P: AsRef<Path>,
{
    // Construct reader and skip initial 5 lines.
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
        let record: HiPrediction =
            record.map_err(|e| anyhow::anyhow!("problem parsing record: {}", e))?;
        result.push(record);
    }

    Ok(result)
}

/// Facade struct that allows easy acccess to the DECIPHER HI data.
#[derive(Debug, Clone)]
pub struct Data {
    /// The DECIPHER HI data.
    data: Vec<HiPrediction>,
    /// Mapping from HGNC identifier to DECIPHER HI prediction index.
    hgnc_to_data_idx: rustc_hash::FxHashMap<String, usize>,
}

impl Data {
    /// Load from file and construct.
    pub fn load<P>(path: P) -> Result<Self, anyhow::Error>
    where
        P: AsRef<Path>,
    {
        let data = load_file(path)?;
        Ok(Self::new(data))
    }

    /// Create a new `Data` object.
    ///
    /// # Arguments
    ///
    /// * `data` - The DECIPHER HI data.
    ///
    /// # Returns
    ///
    /// A new `Data` object.
    pub fn new(data: Vec<HiPrediction>) -> Self {
        let hgnc_to_data_idx = data
            .iter()
            .enumerate()
            .map(|(idx, hi_prediction)| (hi_prediction.hgnc_id.clone(), idx))
            .collect::<rustc_hash::FxHashMap<_, _>>();
        Self {
            data,
            hgnc_to_data_idx,
        }
    }

    /// Get the DECIPHER HI prediction for the given HGNC identifier.
    ///
    /// # Arguments
    ///
    /// * `hgnc_id` - HGNC identifier.
    ///
    /// # Returns
    ///
    /// The DECIPHER HI prediction for the given HGNC identifier, if any.
    pub fn by_hgnc_id(&self, hgnc_id: &str) -> Option<&HiPrediction> {
        self.hgnc_to_data_idx
            .get(hgnc_id)
            .map(|idx| &self.data[*idx])
    }
}

#[cfg(test)]
mod test {
    #[test]
    fn test_load_file() -> Result<(), anyhow::Error> {
        let genes = super::load_file("tests/data/strucvars/decipher_hi_prediction.tsv")?;
        assert_eq!(genes.len(), 17995);

        Ok(())
    }

    #[test]
    fn data_load() -> Result<(), anyhow::Error> {
        let data = super::Data::load("tests/data/strucvars/decipher_hi_prediction.tsv")?;
        insta::assert_yaml_snapshot!(data.by_hgnc_id("HGNC:1100"));

        Ok(())
    }
}
