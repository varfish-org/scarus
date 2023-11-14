//! Implementation of evaluation of copy number loss section 3.

use super::result::Section;
use crate::strucvars::eval::del::result::{L3Count, L1, L1A, L3};

/// Evaluation of deletions, loss of copy number.
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

    /// Perform the evaluation of copy number loss Section 3 and all subsection.
    ///
    /// # Arguments
    ///
    /// * `l1` - Result of section 1.
    ///
    /// # Returns
    ///
    /// The evaluation result.
    ///
    /// # Errors
    ///
    /// If anything goes wrong, it returns a generic `anyhow::Error`.
    pub fn evaluate(&self, l1: &Section) -> Result<Section, anyhow::Error> {
        tracing::debug!("re-evaluating Section 1 results for Section 3 for {:?}", l1);
        let result = match l1 {
            Section::L1(L1::L1A(L1A { genes })) => {
                if genes.len() < 25 {
                    Section::L3(L3::L3A(L3Count {
                        suggested_score: 0.0,
                        num_genes: genes.len(),
                        genes: genes.clone(),
                    }))
                } else if genes.len() < 35 {
                    Section::L3(L3::L3A(L3Count {
                        suggested_score: 0.45,
                        num_genes: genes.len(),
                        genes: genes.clone(),
                    }))
                } else {
                    Section::L3(L3::L3A(L3Count {
                        suggested_score: 0.9,
                        num_genes: genes.len(),
                        genes: genes.clone(),
                    }))
                }
            }
            Section::L1(L1::L1B(_)) => Section::L3(L3::L3A(L3Count {
                suggested_score: 0.0,
                num_genes: 0,
                genes: Vec::new(),
            })),
            _ => anyhow::bail!("unexpected Section 1 result: {:?}", l1),
        };

        tracing::debug!("result = {:?}", &result);

        Ok(result)
    }
}

#[cfg(test)]
mod test {
    use crate::strucvars::data::hgnc::GeneIdInfo;
    use crate::strucvars::eval::del::result::{GeneOverlap, L3Count, Section, L1, L1A, L1B, L3};

    // #[tracing_test::traced_test]
    #[rstest::rstest]
    #[case(0, 0.0)]
    #[case(1, 0.0)]
    #[case(24, 0.0)]
    #[case(25, 0.45)]
    #[case(34, 0.45)]
    #[case(35, 0.9)]
    fn test_evaluate_l1a(
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

        let result = evaluator.evaluate(&Section::L1(L1::L1A(L1A { genes })))?;
        match result {
            Section::L3(L3::L3A(L3Count {
                suggested_score,
                num_genes,
                ..
            }))
            | Section::L3(L3::L3B(L3Count {
                suggested_score,
                num_genes,
                ..
            }))
            | Section::L3(L3::L3C(L3Count {
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

    // #[tracing_test::traced_test]
    #[rstest::rstest]
    fn test_evaluate_l1b() -> Result<(), anyhow::Error> {
        let evaluator = super::Evaluator::new();
        insta::assert_yaml_snapshot!(&evaluator.evaluate(&Section::L1(L1::L1B(L1B {})))?);

        Ok(())
    }
}
