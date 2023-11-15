//! Common functions and data structures for evaluation.

use crate::strucvars::data::hgnc::GeneIdInfo;

/// Scores of the results in individual categories.
pub trait SuggestedScore {
    /// Suggested score for the category.
    fn suggested_score(&self) -> f32;
}

/// Score range for a seciton.
pub trait ScoreRange {
    /// Minimal score for the category.
    fn min_score(&self) -> f32;

    /// Maximum score for the category.
    fn max_score(&self) -> f32;
}

/// Per-gene transcript overlaps as part of `L1A`.
#[derive(Debug, Clone, PartialEq, serde::Deserialize, serde::Serialize)]
pub struct GeneOverlap {
    /// Gene identifiers.
    pub gene: GeneIdInfo,
    /// Transcript identifiers of this gene.
    pub tx_ids: Vec<String>,
}

impl GeneOverlap {
    /// Create a new `GeneOverlap`.
    ///
    /// # Arguments
    ///
    /// * `gene` - Gene identifier.
    /// * `tx_ids` - Transcript identifiers of this gene.
    ///
    /// # Returns
    ///
    /// A new `GeneOverlap`.
    pub fn new(gene: GeneIdInfo, tx_ids: Vec<String>) -> Self {
        Self { gene, tx_ids }
    }
}
