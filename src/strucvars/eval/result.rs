//! Common code for evaluation results of CNVs.

use annonars::gnomad_sv::cli::query::Record as GnomadSvRecord;
use annonars::pbs::clinvar::sv::Record as ClinvarSvRecord;
use intervals_general::bound_pair::BoundPair;

use super::{del, dup};

/// Overall result.
#[derive(Debug, Clone, PartialEq, serde::Deserialize, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Evaluation {
    /// Result of assessing DUPs.
    Dup(dup::result::Evaluation),
    /// Result of assessing DELs.
    Del(del::result::Evaluation),
}

/// Overal assessment of a CNV.
#[derive(
    Debug,
    Clone,
    Copy,
    Default,
    PartialEq,
    strum_macros::EnumIter,
    serde::Deserialize,
    serde::Serialize,
)]
#[serde(rename_all = "snake_case")]
pub enum ClinicalSignificance {
    /// Benign
    Benign,
    /// Likely benign
    LikelyBenign,
    /// Uncertain_significance
    #[default]
    UncertainSignificance,
    /// Likely_pathogenic
    LikelyPathogenic,
    /// Pathogenic
    Pathogenic,
}

/// Score interval to use.
pub type ScoreRange = intervals_general::interval::Interval<f32>;

impl ClinicalSignificance {
    /// Return the scores,
    pub fn range(self) -> ScoreRange {
        match self {
            ClinicalSignificance::Benign => ScoreRange::UnboundedClosedRight { right: 0.99 },
            ClinicalSignificance::LikelyBenign => ScoreRange::LeftHalfOpen {
                bound_pair: BoundPair::new(-0.99, -0.90).expect("invalid interval?"),
            },
            ClinicalSignificance::UncertainSignificance => ScoreRange::Open {
                bound_pair: BoundPair::new(-0.90, 0.90).expect("invalid interval?"),
            },
            ClinicalSignificance::LikelyPathogenic => ScoreRange::LeftHalfOpen {
                bound_pair: BoundPair::new(0.90, 0.99).expect("invalid interval?"),
            },
            ClinicalSignificance::Pathogenic => ScoreRange::UnboundedClosedLeft { left: 0.99 },
        }
    }
}

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
#[serde(rename_all = "snake_case")]
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

/// Information about one overlapping ClinVar records.
#[derive(Debug, Clone, PartialEq, serde::Deserialize, serde::Serialize)]
pub struct ClinvarSvOverlap {
    /// Reciprocal overlap.
    pub overlap: f32,
    /// Overlapping ClinVar record.
    pub record: ClinvarSvRecord,
}

/// Information about one overlapping gnomAD records.
#[derive(Debug, Clone, PartialEq, serde::Deserialize, serde::Serialize)]
pub struct GnomadSvOverlap {
    /// Reciprocal overlap.
    pub overlap: f32,
    /// Overlapping ClinVar record.
    pub record: GnomadSvRecord,
}
