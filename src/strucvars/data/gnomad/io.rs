//! I/O code for gnomAD constraints data.

use std::{io::BufReader, path::Path};

use serde::{Deserialize, Deserializer, Serialize, Serializer};

/// Deserialize `Option::None` as `"NA"`.
///
/// cf. https://stackoverflow.com/a/56384732/84349
fn deserialize_option_na<'de, D, T: Deserialize<'de>>(
    deserializer: D,
) -> Result<Option<T>, D::Error>
where
    D: Deserializer<'de>,
{
    // We define a local enum type inside of the function because it is untagged, serde will
    // deserialize as the first variant that it can.
    #[derive(Deserialize)]
    #[serde(untagged)]
    enum MaybeNA<U> {
        // If it can be parsed as Option<T>, it will be..
        Value(Option<U>),
        // ... otherwise try parsing as a string.
        NAString(String),
    }

    // Deserialize into local enum.
    let value: MaybeNA<T> = Deserialize::deserialize(deserializer)?;
    match value {
        // If parsed as T or None, return that.
        MaybeNA::Value(value) => Ok(value),

        // Otherwise, if value is string an "n/a", return None (and fail if it is any other
        // string)
        MaybeNA::NAString(string) => {
            if string == "NA" {
                Ok(None)
            } else {
                Err(serde::de::Error::custom("Unexpected string"))
            }
        }
    }
}

/// Serialize `Option::None` as `"NA"`.
fn serialize_option_na<S, T>(x: &Option<T>, s: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
    T: Serialize,
{
    match x {
        Some(x) => s.serialize_some(x),
        None => s.serialize_str("NA"),
    }
}

/// A record from the gnomAD constraints database.
#[serde_with::skip_serializing_none]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Record {
    /// The Ensembl gene ID.
    pub ensembl_gene_id: String,
    /// The NCBI gene ID.
    pub entrez_id: String,
    /// The HGNC gene symbol.
    pub gene_symbol: String,
    /// The expected number of loss-of-function variants.
    #[serde(
        serialize_with = "serialize_option_na",
        deserialize_with = "deserialize_option_na"
    )]
    pub exp_lof: Option<f64>,
    /// The expected number of missense variants.
    #[serde(
        serialize_with = "serialize_option_na",
        deserialize_with = "deserialize_option_na"
    )]
    pub exp_mis: Option<f64>,
    /// The expected number of synonymous variants.
    #[serde(
        serialize_with = "serialize_option_na",
        deserialize_with = "deserialize_option_na"
    )]
    pub exp_syn: Option<f64>,
    /// The missense-related Z-score.
    #[serde(
        serialize_with = "serialize_option_na",
        deserialize_with = "deserialize_option_na"
    )]
    pub mis_z: Option<f64>,
    /// The observed number of loss-of-function variants.
    #[serde(
        serialize_with = "serialize_option_na",
        deserialize_with = "deserialize_option_na"
    )]
    pub obs_lof: Option<u32>,
    /// The observed number of missense variants.
    #[serde(
        serialize_with = "serialize_option_na",
        deserialize_with = "deserialize_option_na"
    )]
    pub obs_mis: Option<u32>,
    /// The observed number of synonymous variants.
    #[serde(
        serialize_with = "serialize_option_na",
        deserialize_with = "deserialize_option_na"
    )]
    pub obs_syn: Option<u32>,
    /// The loss-of-function observed/expected ratio.
    #[serde(
        serialize_with = "serialize_option_na",
        deserialize_with = "deserialize_option_na"
    )]
    pub oe_lof: Option<f64>,
    /// The lower bound of the loss-of-function observed/expected ratio.
    #[serde(
        serialize_with = "serialize_option_na",
        deserialize_with = "deserialize_option_na"
    )]
    pub oe_lof_lower: Option<f64>,
    /// The upper bound of the loss-of-function observed/expected ratio.
    #[serde(
        serialize_with = "serialize_option_na",
        deserialize_with = "deserialize_option_na"
    )]
    pub oe_lof_upper: Option<f64>,
    /// The missense observed/expected ratio.
    #[serde(
        serialize_with = "serialize_option_na",
        deserialize_with = "deserialize_option_na"
    )]
    pub oe_mis: Option<f64>,
    /// The lower bound of the missense observed/expected ratio.
    #[serde(
        serialize_with = "serialize_option_na",
        deserialize_with = "deserialize_option_na"
    )]
    pub oe_mis_lower: Option<f64>,
    /// The upper bound of the missense observed/expected ratio.
    #[serde(
        serialize_with = "serialize_option_na",
        deserialize_with = "deserialize_option_na"
    )]
    pub oe_mis_upper: Option<f64>,
    /// The synonymous observed/expected ratio.
    #[serde(
        serialize_with = "serialize_option_na",
        deserialize_with = "deserialize_option_na"
    )]
    pub oe_syn: Option<f64>,
    /// The lower bound of the synonymous observed/expected ratio.
    #[serde(
        serialize_with = "serialize_option_na",
        deserialize_with = "deserialize_option_na"
    )]
    pub oe_syn_lower: Option<f64>,
    /// The upper bound of the synonymous observed/expected ratio.
    #[serde(
        serialize_with = "serialize_option_na",
        deserialize_with = "deserialize_option_na"
    )]
    pub oe_syn_upper: Option<f64>,
    /// The probability of loss-of-function intolerance (pLI score).
    #[serde(
        alias = "pLI",
        serialize_with = "serialize_option_na",
        deserialize_with = "deserialize_option_na"
    )]
    pub pli: Option<f64>,
    /// The synonymous-related Z-score.
    #[serde(
        serialize_with = "serialize_option_na",
        deserialize_with = "deserialize_option_na"
    )]
    pub syn_z: Option<f64>,
    /// The probability of loss-of-function intolerance (pLI score) from ExAC.
    #[serde(
        alias = "exac_pLI",
        serialize_with = "serialize_option_na",
        deserialize_with = "deserialize_option_na"
    )]
    pub exac_pli: Option<f64>,
    /// The observed number of loss-of-function variants from ExAC.
    #[serde(
        serialize_with = "serialize_option_na",
        deserialize_with = "deserialize_option_na"
    )]
    pub exac_obs_lof: Option<f64>,
    /// The expected number of loss-of-function variants from ExAC.
    #[serde(
        serialize_with = "serialize_option_na",
        deserialize_with = "deserialize_option_na"
    )]
    pub exac_exp_lof: Option<f64>,
    /// The loss-of-function observed/expected ratio from ExAC.
    #[serde(
        serialize_with = "serialize_option_na",
        deserialize_with = "deserialize_option_na"
    )]
    pub exac_oe_lof: Option<f64>,
}

/// Load gnomAD constraints file.
///
/// # Arguments
///
/// * `path` - Path to gnomAD constraints file.
///
/// # Returns
///
/// gnomAD constraints.
///
/// # Errors
///
/// If anything goes wrong, it returns a generic `anyhow::Error`.
pub fn load_file<P>(path: P) -> Result<Vec<Record>, anyhow::Error>
where
    P: AsRef<Path>,
{
    // Construct reader and skip initial 5 lines.
    tracing::debug!("opening file: {:?}", path.as_ref());
    let reader = std::fs::File::open(path)
        .map_err(|e| anyhow::anyhow!("problem opening file: {}", e))
        .map(BufReader::new)?;

    tracing::debug!("reading records...");
    let mut csv_reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .flexible(true)
        .from_reader(reader);
    let mut result = Vec::new();
    for record in csv_reader.deserialize() {
        let record: Record =
            record.map_err(|e| anyhow::anyhow!("problem parsing record: {}", e))?;
        result.push(record);
    }
    tracing::debug!("read a total of {} records", result.len());

    Ok(result)
}

#[cfg(test)]
mod test {
    #[test]
    fn test_load_file() -> Result<(), anyhow::Error> {
        let genes = super::load_file("tests/data/strucvars/gnomad_constraints.tsv")?;

        assert_eq!(genes.len(), 18481);
        insta::assert_yaml_snapshot!(&genes[0..5]);

        Ok(())
    }
}
