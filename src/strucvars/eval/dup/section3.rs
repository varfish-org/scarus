//! Implementation of evaluation of copy number gain section 3.

use super::result::Section;
use crate::strucvars::eval::dup::result::{G3Count, G1, G1A, G3};

/// Evaluation of duplications, gain of copy number.
///
/// This is mainly used to encapsulate the functionality.  Creating new such
/// objects is very straightforward and cheap.
#[derive(Default)]
pub struct Evaluator {}

impl Evaluator {
    /// Create a new `Evaluator`.
    pub fn new() -> Self {
        Self {}
    }

    /// Perform the evaluation of copy number gain Section 3 and all subsection.
    ///
    /// # Arguments
    ///
    /// * `g1` - Result of section 1.
    ///
    /// # Returns
    ///
    /// The evaluation result.
    ///
    /// # Errors
    ///
    /// If anything goes wrong, it returns a generic `anyhow::Error`.
    pub fn evaluate(&self, g1: &Section) -> Result<Section, anyhow::Error> {
        tracing::debug!("re-evaluating Section 1 results for Section 3 for {:?}", g1);
        let result = match g1 {
            Section::G1(G1::G1A(G1A { genes })) => {
                if genes.len() < 35 {
                    Section::G3(G3::G3A(G3Count {
                        suggested_score: 0.0,
                        num_genes: genes.len(),
                        genes: genes.clone(),
                    }))
                } else if genes.len() < 50 {
                    Section::G3(G3::G3A(G3Count {
                        suggested_score: 0.45,
                        num_genes: genes.len(),
                        genes: genes.clone(),
                    }))
                } else {
                    Section::G3(G3::G3A(G3Count {
                        suggested_score: 0.9,
                        num_genes: genes.len(),
                        genes: genes.clone(),
                    }))
                }
            }
            Section::G1(G1::G1B(_)) => Section::G3(G3::G3A(G3Count {
                suggested_score: 0.0,
                num_genes: 0,
                genes: Vec::new(),
            })),
            _ => anyhow::bail!("unexpected Section 1 result: {:?}", g1),
        };

        tracing::debug!("result = {:?}", &result);

        Ok(result)
    }
}

#[cfg(test)]
mod test {
    use crate::strucvars::data::hgnc::GeneIdInfo;
    use crate::strucvars::eval::common::GeneOverlap;
    use crate::strucvars::eval::dup::result::{G3Count, Section, G1, G1A, G1B, G3};

    #[tracing_test::traced_test]
    #[rstest::rstest]
    #[case(0, 0.0)]
    #[case(1, 0.0)]
    #[case(34, 0.0)]
    #[case(35, 0.45)]
    #[case(49, 0.45)]
    #[case(50, 0.9)]
    fn test_evaluate_g1a(
        #[case] n_genes: usize,
        #[case] expected_score: f32,
    ) -> Result<(), anyhow::Error> {
        mehari::common::set_snapshot_suffix!("{}", n_genes);

        let genes: Vec<GeneOverlap> = (0..n_genes)
            .map(|_| GeneOverlap {
                gene: GeneIdInfo::default(),
                tx_ids: vec![String::from("fake-tx-id")],
            })
            .collect();

        let evaluator = super::Evaluator::new();

        let result = evaluator.evaluate(&Section::G1(G1::G1A(G1A { genes })))?;
        match result {
            Section::G3(G3::G3A(G3Count {
                suggested_score,
                num_genes,
                ..
            }))
            | Section::G3(G3::G3B(G3Count {
                suggested_score,
                num_genes,
                ..
            }))
            | Section::G3(G3::G3C(G3Count {
                suggested_score,
                num_genes,
                ..
            })) => {
                assert_eq!(expected_score, suggested_score);
                assert_eq!(n_genes, num_genes);
            }
            _ => panic!("invalid section {:?}", &result),
        }

        Ok(())
    }

    #[tracing_test::traced_test]
    #[rstest::rstest]
    fn test_evaluate_g1b() -> Result<(), anyhow::Error> {
        let evaluator = super::Evaluator::new();
        insta::assert_yaml_snapshot!(&evaluator.evaluate(&Section::G1(G1::G1B(G1B {})))?);

        Ok(())
    }
}
