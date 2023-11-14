//! Common functions and data strucdtures for evaluation.

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
