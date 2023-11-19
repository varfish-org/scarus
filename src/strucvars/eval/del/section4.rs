//! Implementation of evaluation of copy number loss section 4.
//!
//! Note that we can only determine 4O reliably in an automated fashion.  We support
//! reporting overlapping variants from ClinVar but leave them as "dangle" because a
//! human must evaluate the phenotype.

use super::result::Section;
use crate::strucvars::ds::StructuralVariant;

/// Evaluation of deletions, loss of copy number.
///
/// This is mainly used to encapsulate the functionality.  Creating new such
/// objects is very straightforward and cheap.
pub struct Evaluator<'a> {
    /// The parent evaluator.
    #[allow(dead_code)] // TODO
    parent: &'a super::super::Evaluator,
}

impl<'a> Evaluator<'a> {
    /// Create a new `Evaluator`.
    pub fn with_parent(parent: &'a super::super::Evaluator) -> Self {
        Self { parent }
    }

    /// Perform the evaluation of copy number loss Section 3 and all subsection.
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
    pub fn evaluate(&self, strucvar: &StructuralVariant) -> Result<Vec<Section>, anyhow::Error> {
        let _ = strucvar;
        tracing::warn!("Section 4 evaluation not implemented yet");
        Ok(Vec::new())
    }
}

#[cfg(test)]
mod test {
    use crate::strucvars::ds::StructuralVariant;

    use super::super::super::test::global_evaluator_37;

    // #[tracing_test::traced_test]
    #[rstest::rstest]
    fn test_evaluate(
        global_evaluator_37: &super::super::super::Evaluator,
    ) -> Result<(), anyhow::Error> {
        let evaluator = super::Evaluator::with_parent(&global_evaluator_37);

        // TODO: write test after we have an actual implementation
        assert_eq!(
            evaluator.evaluate(&StructuralVariant::default())?,
            Vec::new()
        );

        Ok(())
    }
}
