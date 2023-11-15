//! Common code for evaluation results of CNVs.

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
