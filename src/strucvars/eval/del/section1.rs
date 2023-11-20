//! Implementation of evaluation of copy number loss section 1.

use super::result::{Section, L1, L1A, L1B};
use crate::strucvars::{ds::StructuralVariant, eval::common::FunctionalElement};

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
        let functional_elements = self
            .parent
            .functional_overlaps(&strucvar.chrom, strucvar.start, strucvar.stop)
            .map_err(|e| {
                anyhow::anyhow!(
                    "issue with overlap computation of {:?} with functional elements: {}",
                    strucvar,
                    e
                )
            })?
            .into_iter()
            .map(FunctionalElement::Refseq)
            .collect::<Vec<_>>();

        if genes.is_empty() && functional_elements.is_empty() {
            Ok(Section::L1(L1::L1B(L1B::default())))
        } else {
            Ok(Section::L1(L1::L1A(L1A {
                genes,
                functional_elements,
            })))
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
    #[case("1", 8_412_464, 8_877_699, "gene-RERE")]
    #[case("11", 125_423_939, 125_424_204, "no-gene")]
    fn evaluate_l1a(
        #[case] chrom: &str,
        #[case] start: u32,
        #[case] stop: u32,
        #[case] label: &str,
        global_evaluator_37: &super::super::super::Evaluator,
    ) -> Result<(), anyhow::Error> {
        mehari::common::set_snapshot_suffix!("{}", label);

        let evaluator = Evaluator::with_parent(global_evaluator_37);
        let strucvar = ds::StructuralVariant {
            chrom: chrom.to_string(),
            start,
            stop,
            svtype: ds::SvType::Del,
            ambiguous_range: None,
        };

        let res = evaluator.evaluate(&strucvar)?;
        eprintln!("{:?}", &res);

        assert!(matches!(res, Section::L1(L1::L1A(_))));
        insta::assert_yaml_snapshot!(res);

        Ok(())
    }

    #[rstest::rstest]
    #[case("2", 1, 1, "empty")]
    fn evaluate_l1b(
        #[case] chrom: &str,
        #[case] start: u32,
        #[case] stop: u32,
        #[case] label: &str,
        global_evaluator_37: &super::super::super::Evaluator,
    ) -> Result<(), anyhow::Error> {
        mehari::common::set_snapshot_suffix!("{}", label);

        let evaluator = Evaluator::with_parent(global_evaluator_37);
        let strucvar = ds::StructuralVariant {
            chrom: chrom.to_string(),
            start,
            stop,
            svtype: ds::SvType::Del,
            ambiguous_range: None,
        };

        let res = evaluator.evaluate(&strucvar)?;

        assert!(matches!(res, Section::L1(L1::L1B(_))));
        insta::assert_yaml_snapshot!(res);

        Ok(())
    }
}
