//! Implementation of evaluation of copy number loss section 1.

use super::result::{Section, L1, L1A, L1B};
use crate::strucvars::ds::StructuralVariant;

/// Evaluation of deletions, loss of copy number.
///
/// This is mainly used to encapsulate the functionality.  Creating new such
/// objects is very straightforward and cheap.
pub struct Evaluator<'a> {
    /// The parent evaluator.
    parent: &'a super::super::Evaluator,
}

impl<'a> Evaluator<'a> {
    /// Create a new `Evaluator`.
    pub fn with_parent(parent: &'a super::super::Evaluator) -> Self {
        Self { parent }
    }

    /// Perform the evaluation of copy number loss Section 1 and all subsection.
    ///
    /// # Arguments
    ///
    /// * `strucvar` - Structural variant to be evaluated.
    ///
    /// # Returns
    ///
    /// Returns the evaluation results of the section
    ///
    /// # Errors
    ///
    /// If anything goes wrong, it returns a generic `anyhow::Error`.
    pub fn evaluate(&self, strucvar: &StructuralVariant) -> Result<Section, anyhow::Error> {
        let genes = self
            .parent
            .overlapping_genes(&strucvar.chrom, strucvar.start, strucvar.stop)
            .map_err(|e| {
                anyhow::anyhow!("issue with overlap computation of {:?}: {}", strucvar, e)
            })?;

        if !genes.is_empty() {
            Ok(Section::L1(L1::L1A(L1A { genes })))
        } else {
            Ok(Section::L1(L1::L1B(L1B::default())))
        }
    }
}

#[cfg(test)]
pub mod test {
    use crate::strucvars::ds;
    use crate::strucvars::eval::del::result::{Section, L1};

    use super::super::super::test::global_evaluator_37;
    use super::Evaluator;

    #[rstest::rstest]
    fn evaluate_l1a(
        global_evaluator_37: super::super::super::Evaluator,
    ) -> Result<(), anyhow::Error> {
        let evaluator = Evaluator::with_parent(&global_evaluator_37);
        let strucvar = ds::StructuralVariant {
            chrom: "1".to_string(),
            start: 8412464,
            stop: 8877699,
            svtype: ds::SvType::Del,
            ambiguous_range: None,
        };

        let res = evaluator.evaluate(&strucvar)?;

        assert!(matches!(res, Section::L1(L1::L1A(_))));
        insta::assert_yaml_snapshot!(res);

        Ok(())
    }

    #[rstest::rstest]
    fn evaluate_l1b(
        global_evaluator_37: super::super::super::Evaluator,
    ) -> Result<(), anyhow::Error> {
        let evaluator = Evaluator::with_parent(&global_evaluator_37);
        let strucvar = ds::StructuralVariant {
            chrom: "22".to_string(),
            start: 1,
            stop: 1,
            svtype: ds::SvType::Del,
            ambiguous_range: None,
        };

        let res = evaluator.evaluate(&strucvar)?;

        assert!(matches!(res, Section::L1(L1::L1B(_))));
        insta::assert_yaml_snapshot!(res);

        Ok(())
    }
}
