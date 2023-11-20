//! Common code for evaluation results of CNVs.

use annonars::pbs::annonars::clinvar::v1::sv::Record as ClinvarSvRecord;

/// Enumeration describing the PVS1 results.
#[derive(
    Debug,
    Clone,
    Copy,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Default,
    serde::Deserialize,
    serde::Serialize,
)]
pub enum Pvs1Result {
    /// PVS1
    #[default]
    Pvs1,
    /// PVS1_Strong
    Pvs1Strong,
    /// PVS1_Moderate
    Pvs1Moderate,
    /// PVS1_Supporting
    Pvs1Supporting,
}

/// Information about one overlapping ClinVar record.
#[derive(Debug, Clone, PartialEq, serde::Deserialize, serde::Serialize)]
pub struct ClinvarSvOverlap {
    /// Reciprocal overlap.
    pub overlap: f32,
    /// Overlapping ClinVar record.
    pub record: ClinvarSvRecord,
}
