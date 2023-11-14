//! Evaluation of duplications, gain of copy number.

pub mod result;
pub mod section1;
pub mod section2;
pub mod section3;
pub mod section4;

use crate::strucvars::ds::StructuralVariant;

/// Evaluation of duplication, gain of copy number.
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
        // Evaluate section 1: "Contains known functionally important elements".
        // let result_s1 = section1::Evaluator::with_parent(self.parent).evaluate(strucvar)?;
        let result_s1 = section1::Evaluator::with_parent(self.parent).evaluate(strucvar)?;
        // Evaluate section 2: "Overlap wih established/predicted haploinsufficient (HI) or
        // established benign genes/genomic regions (skip to seciton 3 if your copy number loss
        // does not overlap these types of genes/regions)".
        let result_s2 = section2::Evaluator::with_parent(self.parent).evaluate(strucvar)?;
        // Evaluate section 3: "Evaluation of gene number".
        let result_s3 = section3::Evaluator::new().evaluate(&result_s1)?;
        // Evaluate section 4: "Detailed evaluation of genomic content using cases from published literature,
        // public databases, and/or internal lab data".
        let result_s4 = section4::Evaluator::with_parent(self.parent).evaluate(strucvar)?;

        let mut result = Vec::default();
        result.push(result_s1);
        result.append(&mut result_s2.clone());
        result.push(result_s3);
        result.append(&mut result_s4.clone());
        Ok(result)
    }
}

#[cfg(test)]
mod test {
    use super::super::test::global_evaluator_37;
    use crate::strucvars::ds::{StructuralVariant, SvType};
    use crate::strucvars::eval::common::SuggestedScore;

    #[tracing_test::traced_test]
    #[rstest::rstest]
    #[case("1", 11_937_319, 11_944_386, -0.6)] // empty region left of MFN2
    #[case("1", 12_098_550, 12_103_898, -0.6)] // empty region right of MFN2
    #[case("1", 12_050_913, 12_054_733, 0.0)] // contains exon 4 of MFN2
    fn test_evaluate(
        #[case] chrom: &str,
        #[case] start: u32,
        #[case] stop: u32,
        #[case] expected_score: f32,
        global_evaluator_37: super::super::Evaluator,
    ) -> Result<(), anyhow::Error> {
        mehari::common::set_snapshot_suffix!("{}-{}-{}", chrom, start, stop);

        let evaluator = super::Evaluator::with_parent(&global_evaluator_37);
        let result = evaluator.evaluate(&StructuralVariant {
            chrom: chrom.into(),
            start,
            stop,
            svtype: SvType::Dup,
            ..Default::default()
        })?;

        assert_eq_float::assert_eq_float!(
            result[0].suggested_score(),
            expected_score,
            f32::EPSILON
        );
        insta::assert_yaml_snapshot!(result);

        Ok(())
    }
}
