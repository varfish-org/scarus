//! Evaluation of duplications, gain of copy number.

pub mod result;
pub mod section1;
pub mod section3;
pub mod section4;

use crate::strucvars::ds::StructuralVariant;

use self::result::{G1, G1A};

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
        // Evaluate section 1: "Contains known functionally important elements".
        // let result_s1 = section1::Evaluator::with_parent(self.parent).evaluate(strucvar)?;
        let result_s1 = result::Section::G1(G1::G1A(G1A {
            genes: Default::default(),
        }));
        // Evaluate section 2: "Overlap wih established/predicted haploinsufficient (HI) or
        // established benign genes/genomic regions (skip to seciton 3 if your copy number loss
        // does not overlap these types of genes/regions)".
        // let result_s2 = section2::Evaluator::with_parent(self.parent).evaluate(strucvar)?;
        // Evaluate section 3: "Evaluation of gene number".
        let result_s3 = section3::Evaluator::new().evaluate(&result_s1)?;
        // Evaluate section 4: "Detailed evaluation of genomic content using cases from published literature,
        // public databases, and/or internal lab data".
        let result_s4 = section4::Evaluator::with_parent(self.parent).evaluate(strucvar)?;

        let mut result = Vec::default();
        result.push(result_s1);
        // result.append(&mut result_s2.clone());
        result.push(result_s3);
        result.append(&mut result_s4.clone());
        Ok(result)
    }
}
