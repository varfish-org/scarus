//! I/O implementation for `ClinGen` dosage sensitivity.

use std::{
    io::{BufRead, BufReader},
    path::Path,
};

use crate::common::Assembly;

/// Helper trait to adjust the region to the chromosome.
pub trait AdjustForAssembly {
    /// Adjust the region to the chromosome.
    fn for_assembly(self, assembly: Assembly) -> Self;
}

/// `ClinGen` gene dosage sensitivity TSV entry.
#[derive(Debug, PartialEq, serde::Deserialize, serde::Serialize)]
pub struct Gene {
    /// Gene symbol.
    #[serde(alias = "#Gene Symbol")]
    pub gene_symbol: String,
    /// NCBI gene ID.
    #[serde(alias = "Gene ID")]
    pub ncbi_gene_id: String,
    /// Genomic location.
    #[serde(alias = "Genomic Location")]
    pub genomic_location: String,
    /// Haploinsufficiency score.
    #[serde(alias = "Haploinsufficiency Score", deserialize_with = "parse_score")]
    pub haploinsufficiency_score: Option<u32>,
    /// Triplosensitivity score.
    #[serde(alias = "Triplosensitivity Score", deserialize_with = "parse_score")]
    pub triplosensitivity_score: Option<u32>,
    /// Haploinsufficiency Disease ID.
    #[serde(alias = "Haploinsufficiency Disease ID")]
    pub haploinsufficiency_disease_id: Option<String>,
    /// Haploinsufficiency Disease ID.
    #[serde(alias = "Triplosensitivity Disease ID")]
    pub triplosensitivity_disease_id: Option<String>,
}

impl TryInto<super::Gene> for Gene {
    type Error = anyhow::Error;

    fn try_into(self) -> Result<super::Gene, Self::Error> {
        Ok(super::Gene {
            gene_symbol: self.gene_symbol,
            ncbi_gene_id: self.ncbi_gene_id,
            genomic_location: self.genomic_location,
            haploinsufficiency_score: super::Score::try_from(self.haploinsufficiency_score)
                .map_err(|e| anyhow::anyhow!("invalid haploinsufficiency score: {}", e))?,
            triplosensitivity_score: super::Score::try_from(self.triplosensitivity_score)
                .map_err(|e| anyhow::anyhow!("invalid haploinsufficiency score: {}", e))?,
            haploinsufficiency_disease_id: self.haploinsufficiency_disease_id,
            triplosensitivity_disease_id: self.triplosensitivity_disease_id,
        })
    }
}

impl AdjustForAssembly for Gene {
    fn for_assembly(self, assembly: Assembly) -> Self {
        Self {
            genomic_location: if assembly == Assembly::Grch37 {
                self.genomic_location.replace("chr", "")
            } else {
                self.genomic_location
            },
            ..self
        }
    }
}

/// `ClinGen` region dosage sensitivity entry.
#[derive(Debug, PartialEq, serde::Deserialize, serde::Serialize)]
pub struct Region {
    /// ISCA ID
    #[serde(alias = "#ISCA ID")]
    pub isca_id: String,
    /// ISCA Region Name
    #[serde(alias = "ISCA Region Name")]
    pub isca_region_name: String,
    /// Genomic location.
    #[serde(alias = "Genomic Location")]
    pub genomic_location: String,
    /// Haploinsufficiency score.
    #[serde(alias = "Haploinsufficiency Score", deserialize_with = "parse_score")]
    pub haploinsufficiency_score: Option<u32>,
    /// Triplosensitivity score.
    #[serde(alias = "Triplosensitivity Score", deserialize_with = "parse_score")]
    pub triplosensitivity_score: Option<u32>,
    /// Haploinsufficiency Disease ID.
    #[serde(alias = "Haploinsufficiency Disease ID")]
    pub haploinsufficiency_disease_id: Option<String>,
    /// Haploinsufficiency Disease ID.
    #[serde(alias = "Triplosensitivity Disease ID")]
    pub triplosensitivity_disease_id: Option<String>,
}

impl TryInto<super::Region> for Region {
    type Error = anyhow::Error;

    fn try_into(self) -> Result<super::Region, Self::Error> {
        Ok(super::Region {
            isca_id: self.isca_id,
            isca_region_name: self.isca_region_name,
            genomic_location: self.genomic_location,
            haploinsufficiency_score: super::Score::try_from(self.haploinsufficiency_score)
                .map_err(|e| anyhow::anyhow!("invalid haploinsufficiency score: {}", e))?,
            triplosensitivity_score: super::Score::try_from(self.triplosensitivity_score)
                .map_err(|e| anyhow::anyhow!("invalid haploinsufficiency score: {}", e))?,
            haploinsufficiency_disease_id: self.haploinsufficiency_disease_id,
            triplosensitivity_disease_id: self.triplosensitivity_disease_id,
        })
    }
}

impl AdjustForAssembly for Region {
    fn for_assembly(self, assembly: Assembly) -> Self {
        Self {
            genomic_location: if assembly == Assembly::Grch37 {
                self.genomic_location.replace("chr", "")
            } else {
                self.genomic_location
            },
            ..self
        }
    }
}

/// Helper for parsing the scores which may have interesting values.
fn parse_score<'de, D>(d: D) -> Result<Option<u32>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    let tmp: String = serde::Deserialize::deserialize(d)?;
    if tmp.is_empty() || tmp == "Not yet evaluated" || tmp == "-1" {
        Ok(None)
    } else {
        Ok(Some(tmp.parse().map_err(serde::de::Error::custom)?))
    }
}

/// Load `ClinGen` regions file.
///
/// # Arguments
///
/// * `path` - Path to `ClinGen` regions file.
/// * `assembly` - Assembly of the regions file, used to determine whether or not to
///   strip the `chr` prefix from the chromosome name.
///
/// # Returns
///
/// `ClinGen` regions.
///
/// # Errors
///
/// If anything goes wrong, it returns a generic `anyhow::Error`.
pub fn load_file<T, P>(path: P, assembly: Assembly) -> Result<Vec<T>, anyhow::Error>
where
    T: serde::de::DeserializeOwned + std::fmt::Debug + AdjustForAssembly,
    P: AsRef<Path>,
{
    // Construct reader and skip initial 5 lines.
    let mut reader = std::fs::File::open(path)
        .map_err(|e| anyhow::anyhow!("problem opening file: {}", e))
        .map(BufReader::new)?;

    {
        let mut buf = String::new();
        for _ in 0..5 {
            reader.read_line(&mut buf)?;
            buf.clear();
        }
    }

    let mut csv_reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .flexible(true)
        .from_reader(reader);
    let mut result = Vec::new();
    for record in csv_reader.deserialize() {
        let record: T = record.map_err(|e| anyhow::anyhow!("problem parsing record: {}", e))?;
        result.push(record.for_assembly(assembly));
    }

    Ok(result)
}

#[cfg(test)]
mod test {
    use crate::common::Assembly;

    #[test]
    fn load_clingen_genes_file_grch37() -> Result<(), anyhow::Error> {
        let genes = super::load_file::<super::Gene, _>(
            "tests/data/strucvars/ClinGen_gene_curation_list_GRCh37.tsv",
            Assembly::Grch37,
        )?;

        assert_eq!(genes.len(), 1518);
        insta::assert_yaml_snapshot!(&genes[0..5]);

        Ok(())
    }

    #[test]
    fn load_clingen_regions_file_grch37() -> Result<(), anyhow::Error> {
        let regions = super::load_file::<super::Region, _>(
            "tests/data/strucvars/ClinGen_region_curation_list_GRCh37.tsv",
            Assembly::Grch37,
        )?;

        assert_eq!(regions.len(), 513);
        insta::assert_yaml_snapshot!(&regions[0..5]);

        Ok(())
    }

    #[test]
    fn load_clingen_regions_file_grch38() -> Result<(), anyhow::Error> {
        let regions = super::load_file::<super::Region, _>(
            "tests/data/strucvars/ClinGen_region_curation_list_GRCh38.tsv",
            Assembly::Grch38,
        )?;

        assert_eq!(regions.len(), 514);
        insta::assert_yaml_snapshot!(&regions[0..5]);

        Ok(())
    }

    #[test]
    fn gene_conversion() -> Result<(), anyhow::Error> {
        let io_gene = super::Gene {
            gene_symbol: "AASS".to_string(),
            ncbi_gene_id: "53947".to_string(),
            genomic_location: "chr7:121713598-121784344".to_string(),
            haploinsufficiency_score: Some(30),
            triplosensitivity_score: None,
            haploinsufficiency_disease_id: Some("MONDO:0009388".to_string()),
            triplosensitivity_disease_id: None,
        };

        let gene: super::super::Gene = io_gene.try_into()?;

        assert_eq!(
            gene,
            super::super::Gene {
                gene_symbol: "AASS".to_string(),
                ncbi_gene_id: "53947".to_string(),
                genomic_location: "chr7:121713598-121784344".to_string(),
                haploinsufficiency_score: super::super::Score::GeneAssociatedWithRecessivePhenotype,
                triplosensitivity_score: super::super::Score::NoEvidenceAvailable,
                haploinsufficiency_disease_id: Some("MONDO:0009388".to_string()),
                triplosensitivity_disease_id: None,
            }
        );

        Ok(())
    }

    #[test]
    fn region_conversion() -> Result<(), anyhow::Error> {
        let io_region = super::Region {
            isca_id: "ISCA-46748".to_string(),
            isca_region_name: "Xq25 region (includes STAG2)".to_string(),
            genomic_location: "chrX:123034319-123236519".to_string(),
            haploinsufficiency_score: Some(3),
            triplosensitivity_score: Some(3),
            haploinsufficiency_disease_id: Some("MONDO:0100038".to_string()),
            triplosensitivity_disease_id: None,
        };

        let region: super::super::Region = io_region.try_into()?;

        assert_eq!(
            region,
            super::super::Region {
                isca_id: "ISCA-46748".to_string(),
                isca_region_name: "Xq25 region (includes STAG2)".to_string(),
                genomic_location: "chrX:123034319-123236519".to_string(),
                haploinsufficiency_score: super::super::Score::SufficientEvidence,
                triplosensitivity_score: super::super::Score::SufficientEvidence,
                haploinsufficiency_disease_id: Some("MONDO:0100038".to_string()),
                triplosensitivity_disease_id: None,
            }
        );

        Ok(())
    }
}
