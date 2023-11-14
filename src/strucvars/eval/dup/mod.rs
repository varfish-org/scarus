//! Evaluation of duplications, gain of copy number.

pub mod result;

use crate::strucvars::ds::StructuralVariant;

/// Evaluation of deletions, loss of copy number.
///
/// This is mainly used to encapsulate the functionality.  Creating new such
/// objects is very straightforward and cheap.
pub struct Evaluator<'a> {
    /// The parent evaluator.
    parent: &'a super::Evaluator,
}

impl<'a> Evaluator<'a> {
    /// Create a new `Evaluator`.
    pub fn with_parent(parent: &'a super::Evaluator) -> Self {
        Self { parent }
    }

    /// Perform the evaluation of the different sections.
    ///
    /// # Arguments
    ///
    /// * `strucvar` - Structural variant to be evaluated.
    ///
    /// # Returns
    ///
    /// The evaluation result.
    ///
    /// # Errors
    ///
    /// If anything goes wrong, it returns a generic `anyhow::Error`.
    pub fn evaluate(
        &self,
        strucvar: &StructuralVariant,
    ) -> Result<Vec<result::Section>, anyhow::Error> {
        let result = Vec::new();
        Ok(result)
    }
}
