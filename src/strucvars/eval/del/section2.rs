//! Implementation of evaluation of copy number loss section 2.

use annonars::pbs::annonars::clinvar::v1::minimal::{ClinicalSignificance, ReviewStatus};
use bio::bio_types::genome::AbstractInterval;
use hgvs::data::interface::{Provider, TxExonsRecord, TxInfoRecord};

use super::result::{GeneHiPrediction, GeneHiPredictorResult, Section, L2H};
use crate::strucvars::{
    data::{
        clingen_dosage::{Gene, Region, Score},
        intervals::{contains, do_overlap, exon_to_interval, Interval},
    },
    ds::StructuralVariant,
    eval::{common::SuggestedScore as _, del::result::L2D3},
    eval::{
        del::result::{L2, L2A, L2B, L2C1, L2C2, L2D1, L2D2, L2D4, L2E, L2F, L2G},
        result::Pvs1Result,
    },
};

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

    /// Perform the evaluation of copy number loss Section 2 and all subsection.
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
        tracing::debug!("evaluating section 2 for {:?}", strucvar);
        let sv_interval: Interval = strucvar.clone().into();

        // Obtain any overlapping ClinGen region and gene dosage info.
        let (clingen_genes, clingen_regions) = self.parent.clingen_overlaps(&sv_interval);

        // Check each overlapping ClinGen region/gene.  If `sv_interval` completely contains
        // a region/gene and the region/gene has a "sufficient evidence" score for HI then we
        // have case 2A.
        if let Some(result) = Self::handle_case_2a(&clingen_genes, &clingen_regions, &sv_interval) {
            tracing::debug!("case 2A fired: {:?}", &result);
            return Ok(result);
        }
        // Otherwise, we need to test case 2B.
        let has_hi_genes = clingen_genes
            .iter()
            .any(|gene| gene.haploinsufficiency_score == Score::SufficientEvidence);
        let has_hi_region = clingen_regions
            .iter()
            .any(|region| region.haploinsufficiency_score == Score::SufficientEvidence);
        let mut result = if has_hi_genes || has_hi_region {
            let result = vec![Section::L2(L2::L2B(L2B::default()))];
            tracing::debug!("case 2B fired: {:?}", &result);
            result
        } else {
            tracing::debug!("negative for case 2B");
            Vec::new()
        };

        // Handle cases 2C-2E.
        if let Some(result_inner) = self.handle_cases_2c_2e(&sv_interval, &clingen_genes)? {
            result.push(result_inner.clone());
            if !matches!(result_inner, Section::L2(L2::L2D1(L2D1 { .. }))) {
                // Return in all but case 2D-1 (2D-1 = only 3' UTR is involved)
                return Ok(result);
            }
        }

        // Check each overlapping ClinGen region for case 2F: if `sv_interval` is
        // completely contained in this region and the region has a "dosage sensitivity
        // unlikely" score.
        let benign_regions = clingen_regions
            .iter()
            .filter(|clingen_region| {
                clingen_region.haploinsufficiency_score == Score::DosageSensitivityUnlikely
            })
            .cloned()
            .collect::<Vec<_>>();
        if let Some(value) = Self::handle_case_2f(&benign_regions, &sv_interval, &result) {
            return Ok(value);
        }
        // If there is any overlap with a benign region then case 2G is true.
        let result = {
            let mut result = result;
            if !benign_regions.is_empty() {
                // Case 2G positive.
                tracing::debug!("case 2G positive");
                result.push(Section::L2(L2::L2G(L2G {
                    suggested_score: 0.0,
                    benign_regions,
                })));
            } else {
                // Case 2G negative.
                tracing::debug!("case 2G negative");
            };
            result
        };

        // Handle Case 2H: Two or more HI predictors suggest AT LEAST ONE gene in the interval is HI.
        let result = if let Some(result_2h) = self.handle_case_2h(&sv_interval)? {
            let mut result = result;
            result.push(result_2h);
            tracing::debug!("case 2G+2H fired: {:?}", &result);
            result
        } else {
            tracing::debug!("case 2G fired: {:?}", &result);
            result
        };
        Ok(result)
    }

    /// Handle case 2A.
    fn handle_case_2a(
        clingen_genes: &[Gene],
        clingen_regions: &[Region],
        sv_interval: &bio::bio_types::genome::Interval,
    ) -> Option<Vec<Section>> {
        let hi_genes = clingen_genes
            .iter()
            .filter(|clingen_gene| {
                let clingen_interval: Interval = (*clingen_gene)
                    .clone()
                    .try_into()
                    .expect("no interval for gene");
                contains(sv_interval, &clingen_interval)
                    && clingen_gene.haploinsufficiency_score == Score::SufficientEvidence
            })
            .cloned()
            .collect::<Vec<_>>();
        let hi_regions = clingen_regions
            .iter()
            .cloned()
            .filter(|clingen_region| {
                let clingen_interval: Interval = clingen_region
                    .clone()
                    .try_into()
                    .expect("no interval for region");
                contains(sv_interval, &clingen_interval)
                    && clingen_region.haploinsufficiency_score == Score::SufficientEvidence
            })
            .collect::<Vec<_>>();
        if !hi_regions.is_empty() || !hi_genes.is_empty() {
            Some(vec![Section::L2(L2::L2A(L2A {
                suggested_score: 1.0,
                hi_genes,
                hi_regions,
            }))])
        } else {
            None
        }
    }

    // Handle cases 2C-2E.
    fn handle_cases_2c_2e(
        &self,
        sv_interval: &bio::bio_types::genome::Interval,
        clingen_genes: &[Gene],
    ) -> Result<Option<Section>, anyhow::Error> {
        // Get HGNC ids of HI genes.
        let hi_hgnc_ids = clingen_genes
            .iter()
            .filter_map(|clingen_gene| {
                if clingen_gene.haploinsufficiency_score == Score::SufficientEvidence {
                    self.parent
                        .gene_id_data
                        .by_ncbi_gene_id(&clingen_gene.ncbi_gene_id)
                        .map(|gene_id_info| gene_id_info.hgnc_id.clone())
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();
        // Obtain all transcripts of established HI genes.
        if let Some(alt_acc) = self.parent.chrom_to_ac.get(sv_interval.contig()) {
            let tx_ids = self
                .parent
                .provider
                .get_tx_for_region(
                    alt_acc,
                    "splign",
                    sv_interval.range().start as i32,
                    sv_interval.range().end as i32,
                )
                .map_err(|e| anyhow::anyhow!("failed to get transcripts for region: {}", e))?
                .into_iter()
                .filter_map(|tx| {
                    // filter to the transcripts of the HI genes
                    self.parent.provider.get_tx(&tx.tx_ac).and_then(|tx| {
                        if hi_hgnc_ids.contains(&tx.gene_id) {
                            Some(tx.id)
                        } else {
                            None
                        }
                    })
                })
                .collect::<Vec<_>>();
            let mut sections = tx_ids
                .into_iter()
                .filter_map(|tx_id| {
                    let mut exons = self
                        .parent
                        .provider
                        .get_tx_exons(&tx_id, alt_acc, "splign")
                        .expect("no exons for transcript?");
                    // Ensure exons are sorted by reference start position.
                    exons.sort_by_key(|exon| exon.alt_start_i);
                    let tx_info = self
                        .parent
                        .provider
                        .get_tx_info(&tx_id, alt_acc, "splign")
                        .expect("no tx info for transcript?");
                    // Now, test cases 2C..2E.
                    if let Some(section) = self.handle_case_2c(&tx_info, &exons, sv_interval) {
                        Some(section)
                    } else if let Some(section) = self
                        .handle_case_2d(&tx_info, &exons, sv_interval)
                        .ok()
                        .flatten()
                    {
                        Some(section)
                    } else {
                        self.handle_case_2e(&tx_info, sv_interval)
                    }
                })
                .collect::<Vec<_>>();

            // Sort by suggested score.
            sections.sort_by(|a, b| {
                let a = a.suggested_score();
                let b = b.suggested_score();
                a.partial_cmp(&b).unwrap_or(std::cmp::Ordering::Equal)
            });

            Ok(sections.pop())
        } else {
            Ok(None)
        }
    }

    /// Handle case 2C: 5' end overlap.
    fn handle_case_2c(
        &self,
        tx_info: &TxInfoRecord,
        exons: &[TxExonsRecord],
        sv_interval: &bio::bio_types::genome::Interval,
    ) -> Option<Section> {
        let most_lhs = exons
            .iter()
            .map(|exon| exon.alt_start_i - 1)
            .min()
            .expect("no exons?");
        let most_lhs =
            bio::bio_types::genome::Locus::new(sv_interval.contig().to_string(), most_lhs as u64);
        let most_rhs = exons
            .iter()
            .map(|exon| exon.alt_end_i)
            .max()
            .expect("no exons?");
        let most_rhs =
            bio::bio_types::genome::Locus::new(sv_interval.contig().to_string(), most_rhs as u64);
        let is_forward = exons.first().expect("no exons?").alt_strand >= 0;

        if is_forward && !sv_interval.contains(most_lhs.clone()) // forward: 5' is left
            // reverse: 5' is right
            || !is_forward && !sv_interval.contains(most_rhs.clone())
        {
            tracing::debug!(
                "case 2C negative: {}, {:?}, lhs={:?}, rhs={:?}",
                is_forward,
                &sv_interval,
                &most_lhs,
                &most_rhs
            );
            // is not case 2C positive
            return None;
        }
        // Case 2C positive.
        tracing::debug!("case 2C fired: {:?} contains {:?}", &sv_interval, &most_lhs);
        assert!(
            (!is_forward || !sv_interval.contains(most_rhs))
                && (is_forward | !sv_interval.contains(most_lhs))
        );
        if let (Some(cds_start_i), Some(cds_end_i)) = (tx_info.cds_start_i, tx_info.cds_end_i) {
            let cds_interval = Interval::new(
                sv_interval.contig().to_string(),
                (cds_start_i as u64)..(cds_end_i as u64),
            );
            let gene = self
                .parent
                .gene_id_data
                .by_hgnc_id(&tx_info.hgnc)
                .expect("no gene id info?")
                .clone();
            if do_overlap(sv_interval, &cds_interval) {
                // Case 2C-1 positive.
                tracing::debug!("case 2C-1: {:?} vs. {:?}", sv_interval, &cds_interval);
                Some(Section::L2(L2::L2C1(L2C1 {
                    suggested_score: 0.9,
                    gene,
                })))
            } else {
                // Case 2C-2 positive.
                tracing::debug!("case 2C-2: {:?} vs. {:?}", sv_interval, &cds_interval);
                Some(Section::L2(L2::L2C2(L2C2 {
                    suggested_score: 0.0,
                    gene,
                })))
            }
        } else {
            tracing::debug!("non-coding transcript {:?}?", &tx_info.tx_ac);
            None
        }
    }

    /// Handle case 2D: 3' end overlap.
    fn handle_case_2d(
        &self,
        tx_info: &TxInfoRecord,
        exons: &[TxExonsRecord],
        sv_interval: &bio::bio_types::genome::Interval,
    ) -> Result<Option<Section>, anyhow::Error> {
        let most_lhs = exons
            .iter()
            .map(|exon| exon.alt_start_i)
            .min()
            .expect("no exons?");
        let most_lhs =
            bio::bio_types::genome::Locus::new(sv_interval.contig().to_string(), most_lhs as u64);
        let most_rhs = exons
            .iter()
            .map(|exon| exon.alt_end_i - 1)
            .max()
            .expect("no exons?");
        let most_rhs =
            bio::bio_types::genome::Locus::new(sv_interval.contig().to_string(), most_rhs as u64);
        let is_forward = exons.first().expect("no exons?").alt_strand >= 0;

        if is_forward && !sv_interval.contains(most_rhs.clone()) // forward: 5' is right
            || !is_forward && !sv_interval.contains(most_lhs.clone())
        {
            // is not case 2D positive
            return Ok(None);
        }
        // Case 2D positive.
        tracing::debug!("case 2D fired: {:?} contains {:?}", &sv_interval, &most_rhs);
        assert!(
            (!is_forward || !sv_interval.contains(most_lhs))
                && (is_forward | !sv_interval.contains(most_rhs))
        );

        if let (Some(cds_start_i), Some(cds_end_i)) = (tx_info.cds_start_i, tx_info.cds_end_i) {
            // Get CDS as interval.
            let cds_interval = Interval::new(
                sv_interval.contig().to_string(),
                (cds_start_i as u64)..(cds_end_i as u64),
            );

            let gene = self
                .parent
                .gene_id_data
                .by_hgnc_id(&tx_info.hgnc)
                .expect("no gene id info?")
                .clone();
            if !do_overlap(sv_interval, &cds_interval) {
                // Case 2D-1 positive: only 3' UTR is involved.
                tracing::debug!("case 2D-1 fired: {:?} vs. {:?}", sv_interval, &cds_interval);
                return Ok(Some(Section::L2(L2::L2D1(L2D1 { gene }))));
            }

            // Get all exons that are affected and in CDS.
            let affected_cds_exons = exons
                .iter()
                .filter(|exon| {
                    let exon = exon_to_interval(sv_interval.contig().to_string(), exon);
                    do_overlap(&exon, sv_interval) && do_overlap(&exon, &cds_interval)
                })
                .collect::<Vec<_>>();

            assert!(!affected_cds_exons.is_empty());
            if affected_cds_exons.len() == 1 {
                let clinvar_records = self
                    .parent
                    .query_clinvar_range(
                        sv_interval.contig(),
                        affected_cds_exons[0].alt_start_i as u32 + 1,
                        affected_cds_exons[0].alt_end_i as u32,
                    )
                    .map_err(|e| anyhow::anyhow!("failed to query ClinVar database: {}", e))?;
                // We count VCV records as established if they have two stars.
                let mut established_pathogenic_vcvs = clinvar_records
                    .into_iter()
                    .filter_map(|clinvar_record| {
                        let clinsig = clinvar_record.reference_assertions[0].clinical_significance;
                        let status = clinvar_record.reference_assertions[0].review_status;
                        if status
                            <= ReviewStatus::CriteriaProvidedMultipleSubmittersNoConflicts as i32
                            && clinsig <= ClinicalSignificance::LikelyPathogenic as i32
                        {
                            Some(clinvar_record.vcv)
                        } else {
                            None
                        }
                    })
                    .collect::<Vec<_>>();
                established_pathogenic_vcvs.sort();
                if established_pathogenic_vcvs.is_empty() {
                    tracing::debug!("case 2D-3 fired: {:?}", affected_cds_exons);
                    Ok(Some(Section::L2(L2::L2D3(L2D3 {
                        suggested_score: 0.3,
                        gene,
                    }))))
                } else {
                    tracing::debug!("case 2D-2 fired: {:?}", affected_cds_exons);
                    Ok(Some(Section::L2(L2::L2D2(L2D2 {
                        suggested_score: 0.9,
                        gene,
                        clinvar_ids: established_pathogenic_vcvs,
                    }))))
                }
            } else {
                // Case 2D-4 positive: other exons in addition to last one involved.  Nonsense-
                // mediated decay is expected to occur.  Note that we interpret the second
                // sentence as an instruction here and not a second condition.
                tracing::debug!("case 2D-4 fired: {:?}", affected_cds_exons);
                Ok(Some(Section::L2(L2::L2D4(L2D4 {
                    suggested_score: 0.90,
                    gene,
                    exon_nos: affected_cds_exons
                        .iter()
                        .map(|exon| exon.ord as u32 + 1)
                        .collect::<Vec<_>>(),
                }))))
            }
        } else {
            tracing::debug!("non-coding transcript {:?}?", &tx_info.tx_ac);
            Ok(None)
        }
    }

    /// Handle case 2E.
    fn handle_case_2e(
        &self,
        tx_info: &TxInfoRecord,
        sv_interval: &bio::bio_types::genome::Interval,
    ) -> Option<Section> {
        let _ = tx_info;
        let _ = sv_interval;
        // Case 2E positive.
        tracing::debug!("case 2E fired");
        tracing::warn!("PVS1 not yet implemented for case 2E TODO !!! TODO");

        Some(Section::L2(L2::L2E(L2E {
            suggested_score: 0.9,
            pvs1_result: Pvs1Result::Pvs1,
        })))
    }

    /// Handle case 2F: complete contained within benign region.
    fn handle_case_2f(
        benign_regions: &[Region],
        sv_interval: &bio::bio_types::genome::Interval,
        result: &[Section],
    ) -> Option<Vec<Section>> {
        let containing_benign_regions = benign_regions
            .iter()
            .filter(|clingen_region| {
                let clingen_interval: Interval = (*clingen_region)
                    .clone()
                    .try_into()
                    .expect("no interval for region");
                contains(&clingen_interval, sv_interval)
            })
            .map(|region| (*region).clone())
            .collect::<Vec<_>>();
        if !containing_benign_regions.is_empty() {
            let mut result = Vec::from_iter(result.iter().cloned());
            result.push(Section::L2(L2::L2F(L2F {
                suggested_score: -1.0,
                benign_regions: containing_benign_regions,
            })));
            tracing::debug!("case 2F fired: {:?}", &result);
            Some(result)
        } else {
            None
        }
    }

    /// Handle case 2H.
    ///
    /// Evaluate whether any overlapping RefSeq gene has two or more HI predictors
    /// suggesting that it HI.
    ///
    /// We use the criterion from AutoCNV. (1) gnomAD pLI>=0.9 and upper bound of
    /// observed/expected LoF confidence interval <0.35 and (2) DECIPHER HI index
    /// <=10%.
    fn handle_case_2h(
        &self,
        sv_interval: &bio::bio_types::genome::Interval,
    ) -> Result<Option<Section>, anyhow::Error> {
        if let Some(alt_acc) = self.parent.chrom_to_ac.get(sv_interval.contig()) {
            let txs = self
                .parent
                .provider
                .get_tx_for_region(
                    alt_acc,
                    "splign",
                    sv_interval.range().start as i32,
                    sv_interval.range().end as i32,
                )
                .map_err(|e| anyhow::anyhow!("failed to get transcripts for region: {}", e))?;
            let mut hgnc_ids = txs
                .iter()
                .flat_map(|tx| self.parent.provider.get_tx(&tx.tx_ac).map(|tx| tx.gene_id))
                .collect::<Vec<_>>();
            hgnc_ids.sort();
            hgnc_ids.dedup();

            let gene_hi_predictions = hgnc_ids
                .iter()
                .flat_map(|hgnc_id| self.get_gene_hi_predictions(hgnc_id))
                .collect::<Vec<_>>();
            if gene_hi_predictions.is_empty() {
                Ok(None)
            } else {
                Ok(Some(Section::L2(L2::L2H(L2H {
                    suggested_score: 0.15,
                    gene_hi_predictions,
                }))))
            }
        } else {
            Ok(None)
        }
    }

    /// Obtains `GeneHiPrediction`s for the given HGNC ID.
    ///
    /// And checks whether the gene is predicted to be haploinsufficient by two or more
    /// predictors using the criteria outlined in documentation of `Self::handle_case_2h()`.
    fn get_gene_hi_predictions(&self, hgnc_id: &str) -> Option<GeneHiPrediction> {
        let hi_record = self.parent.decipher_hi_data.by_hgnc_id(hgnc_id);
        let gnomad_record = self.parent.gnomad_constraint_data.by_hgnc_id(hgnc_id);
        tracing::info!("HI = {:?}, gnomad = {:?}", &hi_record, &gnomad_record);
        if let (Some(hi_record), Some(gnomad_record)) = (hi_record, gnomad_record) {
            if let (Some(pli), Some(oe_lof_upper)) = (gnomad_record.pli, gnomad_record.oe_lof_upper)
            {
                if hi_record.hi_index <= 10.0 && pli >= 0.9 && oe_lof_upper < 0.35 {
                    Some(GeneHiPrediction {
                        gene: self
                            .parent
                            .gene_id_data
                            .by_hgnc_id(hgnc_id)
                            .expect("could not resolve HGNC ID")
                            .clone(),
                        results: vec![
                            GeneHiPredictorResult::DecipherHiIndex {
                                hi_index: hi_record.hi_index,
                            },
                            GeneHiPredictorResult::GnomadPli {
                                pli_score: pli,
                                oe_lof_upper,
                            },
                        ],
                    })
                } else {
                    None
                }
            } else {
                None
            }
        } else {
            None
        }
    }
}

#[cfg(test)]
mod test {
    use hgvs::data::interface::Provider;

    use crate::strucvars::ds;

    use super::super::super::test::global_evaluator_37;

    /// Test the `evaluate` function with the different cases of Section 2.
    ///
    /// We make one test with `COL3A1` gene on forward strand and `APOB`
    /// on the reverse strand.  We use `COL16A1` as a gene that is unlikely HI
    /// and `REV3L` as one that is only predicted as HI.  `AVCRL2` has pathogenic
    /// variants in the last coding exon.
    ///
    /// We only make checks for `APOB` for cases 2C..2E where the strand is
    /// considered and important for UTR/CDS distinction.
    #[tracing_test::traced_test]
    #[rstest::rstest]
    // Case 2A: with region ISCA-46303
    #[case("17", 67_892_996, 69_792_434, "ISCA-46303", "2A-pos")] // full
    #[case("17", 67_892_996, 69_792_433, "ISCA-46303", "2A-neg-1")] // partial left
    #[case("17", 67_892_997, 69_792_434, "ISCA-46303", "2A-neg-2")] // partial right
    #[case("17", 67_892_997, 69_792_433, "ISCA-46303", "2A-neg-3")]
    // contained
    // Case 2A: with gene COL3A1 (negative cases have no overlap so other cases don't trigger)
    #[case("2", 189_839_099, 189_877_472, "COL3A1", "2A-pos")] // full
    #[case("2", 189_877_473, 189_878_473, "COL3A1", "2A-neg-1")] // left
    #[case("2", 189_838_099, 189_839_098, "COL3A1", "2A-neg-2")] // right
    // Case 2B negative (positive is handled as part of 2C..2E)
    #[case("1", 32_117_848, 32_169_768, "COL16A1", "2B-neg")] // full overlap
    // Case 2C-1: 5' end overlap of HI gene, with coding sequence
    #[case("2", 189_839_029, 189_839_743, "COL3A1", "2C-1")]
    #[case("2", 21_266_811, 21_266_952, "APOB", "2C-1")]
    // Case 2C-2: 5' end overlap of HI gene, only 5' UTR
    #[case("2", 189_839_029, 189_839_210, "COL3A1", "2C-2")]
    #[case("2", 21_266_820, 21_266_952, "APOB", "2C-2")]
    // Case 2D-1: 3' end overlap of HI gene, only 3' UTR
    #[case("2", 189_876_511, 189_877_563, "COL3A1", "2D-1")]
    #[case("2", 21_224_301, 21_224_585, "APOB", "2D-1")]
    // Case 2D-2: 3' end overlap of HI gene, last exon only, pathogenic vars
    #[case("12", 52_314_363, 52_317_600, "AVCRL2", "2D-2")]
    // Case 2D-3: 3' end overlap of HI gene, last exon only, NO pathogenic vars
    #[case("2", 189_875_865, 189_877_563, "COL3A1", "2D-3")]
    #[case("2", 21_224_301, 21_226_760, "APOB", "2D-3")]
    // Case 2D-4: 3' end overlap of HI gene, more than last exon
    #[case("2", 189_875_551, 189_877_563, "COL3A1", "2D-4")]
    #[case("2", 21_224_301, 21_227_294, "APOB", "2D-4")]
    // Case 2E: intragenic variant, both breakpoints in same gene
    // TODO: need more tests once VPS1 has been implemented
    #[case("2", 189_839_528, 189_876_186, "COL3A1", "2E")]
    #[case("2", 21_226_513, 21_266_591, "APOB", "2E")]
    // Case 2F: completely contained in benign region
    #[case("1", 12_989_199, 12_998_900, "ISCA-46311", "2F-pos")] // contained
    #[case("1", 12_998_901, 12_998_901, "ISCA-46311", "2F-neg-1")] // right of
    #[case("1", 12_989_198, 12_989_198, "ISCA-46311", "2F-neg-2")] // left of
    // Case 2G: overlapping benign region, but additional material
    #[case("1", 12_989_199, 12_998_901, "ISCA-46311", "2G-pos-1")] // additional right
    #[case("1", 12_989_198, 12_998_900, "ISCA-46311", "2G-pos-2")]
    // additional left
    // Case 2H: two or more HI predictors suggest AT LEAST ONE gene in the interval is HI
    #[case("1", 12_040_238, 12_073_572, "MFN2", "2H-neg")] // MFN2 not DECIPHER HI
    #[case("6", 111_620_236, 111_804_914, "REV3L", "2H-pos-1")] // REV3L, contained
    #[case("6", 111_620_236, 111_620_236, "REV3L", "2H-pos-2")] // REV3L, 1bp left
    #[case("6", 111_804_914, 111_804_914, "REV3L", "2H-pos-3")] // REV3L, 1bp right
    #[case("6", 111_620_235, 111_620_235, "REV3L", "2H-neg-1")] // REV3L, left of
    #[case("6", 111_804_915, 111_804_915, "REV3L", "2H-neg-2")] // REV3L, right of
    fn evaluate(
        #[case] chrom: &str,
        #[case] start: u64,
        #[case] stop: u64,
        #[case] gene: &str,
        #[case] label: &str,
        global_evaluator_37: &super::super::super::Evaluator,
    ) -> Result<(), anyhow::Error> {
        mehari::common::set_snapshot_suffix!("{}-{}", gene, label);

        let evaluator = super::Evaluator::with_parent(&global_evaluator_37);
        let strucvar = ds::StructuralVariant {
            chrom: chrom.into(),
            start: start as u32,
            stop: stop as u32,
            svtype: ds::SvType::Del,
            ambiguous_range: None,
        };

        let res = evaluator.evaluate(&strucvar)?;
        insta::assert_yaml_snapshot!(res);

        Ok(())
    }

    /// Test inernal working of `handle_case_2a` (complete overlap).
    #[tracing_test::traced_test]
    #[rstest::rstest]
    // Note: same cases as in `evaluate` above -- keep in sync!
    #[case("17", 67_892_996, 69_792_434, "ISCA-46303", "2A-pos")]
    #[case("2", 189_839_099, 189_877_472, "COL3A1", "2A-pos")]
    fn handle_case_2a(
        #[case] chrom: &str,
        #[case] start: u64,
        #[case] stop: u64,
        #[case] hgnc_id: &str,
        #[case] label: &str,
        global_evaluator_37: &super::super::super::Evaluator,
    ) -> Result<(), anyhow::Error> {
        mehari::common::set_snapshot_suffix!("{}-{}", hgnc_id, label);

        let strucvar = ds::StructuralVariant {
            chrom: chrom.into(),
            start: start as u32,
            stop: stop as u32,
            svtype: ds::SvType::Del,
            ambiguous_range: None,
        };
        let sv_interval = strucvar.into();
        let genes = global_evaluator_37
            .clingen_dosage_data
            .gene_by_overlap(&sv_interval);
        let regions = global_evaluator_37
            .clingen_dosage_data
            .region_by_overlap(&sv_interval);

        let res = super::Evaluator::handle_case_2a(&genes, &regions, &sv_interval);
        insta::assert_yaml_snapshot!(res);

        Ok(())
    }

    /// Test internal workings of `handle_cases_2c_2e` (5', 3' overlap, intragenic).
    #[tracing_test::traced_test]
    #[rstest::rstest]
    // Note: same cases as in `evaluate` above -- keep in sync!
    // Case 2C-1: 5' end overlap of HI gene, with coding sequence
    #[case("2", 189_839_029, 189_839_743, "COL3A1", "2C-1")]
    #[case("2", 21_266_811, 21_266_952, "APOB", "2C-1")]
    // Case 2C-2: 5' end overlap of HI gene, only 5' UTR
    #[case("2", 189_839_029, 189_839_210, "COL3A1", "2C-2")]
    #[case("2", 21_266_820, 21_266_952, "APOB", "2C-2")]
    // Case 2D-1: 3' end overlap of HI gene, only 3' UTR
    #[case("2", 189_876_511, 189_877_563, "COL3A1", "2D-1")]
    #[case("2", 21_224_301, 21_224_585, "APOB", "2D-1")]
    // Case 2D-2: 3' end overlap of HI gene, last exon only, pathogenic vars
    #[case("12", 52_314_363, 52_317_600, "AVCRL2", "2D-2")]
    // Case 2D-3: 3' end overlap of HI gene, last exon only, NO pathogenic vars
    #[case("2", 189_875_865, 189_877_563, "COL3A1", "2D-3")]
    #[case("2", 21_224_301, 21_226_760, "APOB", "2D-3")]
    // Case 2D-4: 3' end overlap of HI gene, more than last exon
    #[case("2", 189_875_551, 189_877_563, "COL3A1", "2D-4")]
    #[case("2", 21_224_301, 21_227_294, "APOB", "2D-4")]
    // Case 2E: intragenic variant, both breakpoints in same gene
    // TODO: need more tests once VPS1 has been implemented
    #[case("2", 189_839_528, 189_876_186, "COL3A1", "2E")]
    #[case("2", 21_226_513, 21_266_591, "APOB", "2E")]
    fn handle_cases_2c_2e(
        #[case] chrom: &str,
        #[case] start: u64,
        #[case] stop: u64,
        #[case] hgnc_id: &str,
        #[case] label: &str,
        global_evaluator_37: &super::super::super::Evaluator,
    ) -> Result<(), anyhow::Error> {
        mehari::common::set_snapshot_suffix!("{}-{}", hgnc_id, label);

        let strucvar = ds::StructuralVariant {
            chrom: chrom.into(),
            start: start as u32,
            stop: stop as u32,
            svtype: ds::SvType::Del,
            ambiguous_range: None,
        };
        let sv_interval = strucvar.into();
        let evaluator = super::Evaluator::with_parent(&global_evaluator_37);
        let (genes, _) = global_evaluator_37.clingen_overlaps(&sv_interval);

        let res = evaluator.handle_cases_2c_2e(&sv_interval, &genes)?;
        insta::assert_yaml_snapshot!(res);

        Ok(())
    }

    /// Test internal workings of `handle_case_2c` (5' end overlap from upstream).
    #[tracing_test::traced_test]
    #[rstest::rstest]
    // Note: same cases as in `evaluate` above -- keep in sync!
    // Case 2C-1: 5' end overlap of HI gene, with coding sequence
    #[case("2", 189_839_029, 189_839_743, "COL3A1", "2C-1")]
    #[case("2", 21_266_811, 21_266_952, "APOB", "2C-1")]
    // Case 2C-2: 5' end overlap of HI gene, only 5' UTR
    #[case("2", 189_839_029, 189_839_210, "COL3A1", "2C-2")]
    #[case("2", 21_266_820, 21_266_952, "APOB", "2C-2")]
    fn handle_case_2c(
        #[case] chrom: &str,
        #[case] start: u64,
        #[case] stop: u64,
        #[case] hgnc_id: &str,
        #[case] label: &str,
        global_evaluator_37: &super::super::super::Evaluator,
    ) -> Result<(), anyhow::Error> {
        mehari::common::set_snapshot_suffix!("{}-{}", hgnc_id, label);

        let strucvar = ds::StructuralVariant {
            chrom: chrom.into(),
            start: start as u32,
            stop: stop as u32,
            svtype: ds::SvType::Del,
            ambiguous_range: None,
        };
        let sv_interval = strucvar.into();

        let tx_id = if hgnc_id == "COL3A1" {
            "NM_000090.4"
        } else {
            "NM_000384.3"
        };

        let tx_info = global_evaluator_37
            .provider
            .get_tx_info(tx_id, "NC_000002.11", "splign")?;
        let exons = global_evaluator_37
            .provider
            .get_tx_exons(tx_id, "NC_000002.11", "splign")?;
        let evaluator = super::Evaluator::with_parent(&global_evaluator_37);

        let res = evaluator.handle_case_2c(&tx_info, &exons, &sv_interval);
        insta::assert_yaml_snapshot!(res);

        Ok(())
    }

    /// Test internal workings of `handle_case_2d` (3' end overlap from downstream).
    #[tracing_test::traced_test]
    #[rstest::rstest]
    // Note: same cases as in `evaluate` above -- keep in sync!
    // Case 2D-1: 3' end overlap of HI gene, only 3' UTR
    #[case("2", 189_876_511, 189_877_563, "COL3A1", "2D-1")]
    #[case("2", 21_224_301, 21_224_585, "APOB", "2D-1")]
    // Case 2D-2: 3' end overlap of HI gene, last exon only, pathogenic vars
    // TODO: need to find example with pathogenic variants
    // Case 2D-3: 3' end overlap of HI gene, last exon only, NO pathogenic vars
    #[case("2", 189_875_865, 189_877_563, "COL3A1", "2D-3")]
    #[case("2", 21_224_301, 21_226_760, "APOB", "2D-3")]
    // Case 2D-4: 3' end overlap of HI gene, more than last exon
    #[case("2", 189_875_551, 189_877_563, "COL3A1", "2D-4")]
    #[case("2", 21_224_301, 21_227_294, "APOB", "2D-4")]
    fn handle_case_2d(
        #[case] chrom: &str,
        #[case] start: u64,
        #[case] stop: u64,
        #[case] hgnc_id: &str,
        #[case] label: &str,
        global_evaluator_37: &super::super::super::Evaluator,
    ) -> Result<(), anyhow::Error> {
        mehari::common::set_snapshot_suffix!("{}-{}", hgnc_id, label);

        let strucvar = ds::StructuralVariant {
            chrom: chrom.into(),
            start: start as u32,
            stop: stop as u32,
            svtype: ds::SvType::Del,
            ambiguous_range: None,
        };
        let sv_interval = strucvar.into();

        let tx_id = if hgnc_id == "COL3A1" {
            "NM_000090.4"
        } else {
            "NM_000384.3"
        };

        let tx_info = global_evaluator_37
            .provider
            .get_tx_info(tx_id, "NC_000002.11", "splign")?;
        let exons = global_evaluator_37
            .provider
            .get_tx_exons(tx_id, "NC_000002.11", "splign")?;
        let evaluator = super::Evaluator::with_parent(&global_evaluator_37);

        let res = evaluator.handle_case_2d(&tx_info, &exons, &sv_interval)?;

        insta::assert_yaml_snapshot!(res);

        Ok(())
    }

    /// Test internal workings of `handle_case_2e` (2E positive case).
    #[tracing_test::traced_test]
    #[rstest::rstest]
    // Note: same cases as in `evaluate` above -- keep in sync!
    // Case 2E: intragenic variant, both breakpoints in same gene
    // TODO: need more tests once VPS1 has been implemented
    #[case("2", 189_839_528, 189_876_186, "COL3A1", "2E")]
    #[case("2", 21_226_513, 21_266_591, "APOB", "2E")]
    fn handle_case_2e(
        #[case] chrom: &str,
        #[case] start: u64,
        #[case] stop: u64,
        #[case] hgnc_id: &str,
        #[case] label: &str,
        global_evaluator_37: &super::super::super::Evaluator,
    ) -> Result<(), anyhow::Error> {
        mehari::common::set_snapshot_suffix!("{}-{}", hgnc_id, label);

        let strucvar = ds::StructuralVariant {
            chrom: chrom.into(),
            start: start as u32,
            stop: stop as u32,
            svtype: ds::SvType::Del,
            ambiguous_range: None,
        };
        let sv_interval = strucvar.into();

        let tx_id = if hgnc_id == "COL3A1" {
            "NM_000090.4"
        } else {
            "NM_000384.3"
        };

        let tx_info = global_evaluator_37
            .provider
            .get_tx_info(tx_id, "NC_000002.11", "splign")?;
        let evaluator = super::Evaluator::with_parent(&global_evaluator_37);

        let section = evaluator.handle_case_2e(&tx_info, &sv_interval);
        insta::assert_yaml_snapshot!(section);

        Ok(())
    }

    /// Test internal workings of `handle_case_2f`
    #[tracing_test::traced_test]
    #[rstest::rstest]
    // Note: same cases as in `evaluate` above -- keep in sync!
    // Case 2F: completely contained in benign region
    #[case("1", 12_989_199, 12_998_900, "ISCA-46311", "2F-pos")] // contained
    #[case("1", 12_998_901, 12_998_901, "ISCA-46311", "2F-neg-1")] // right of
    #[case("1", 12_989_198, 12_989_198, "ISCA-46311", "2F-neg-2")] // left of
    fn handle_case_2f(
        #[case] chrom: &str,
        #[case] start: u64,
        #[case] stop: u64,
        #[case] hgnc_id: &str,
        #[case] label: &str,
        global_evaluator_37: &super::super::super::Evaluator,
    ) -> Result<(), anyhow::Error> {
        mehari::common::set_snapshot_suffix!("{}-{}", hgnc_id, label);

        let strucvar = ds::StructuralVariant {
            chrom: chrom.into(),
            start: start as u32,
            stop: stop as u32,
            svtype: ds::SvType::Del,
            ambiguous_range: None,
        };
        let sv_interval = strucvar.into();

        let (gene_overlaps, region_overlaps) = global_evaluator_37.clingen_overlaps(&sv_interval);

        insta::assert_yaml_snapshot!(&gene_overlaps);
        insta::assert_yaml_snapshot!(&region_overlaps);

        let res = super::Evaluator::handle_case_2f(&region_overlaps, &sv_interval, &[]);
        assert_eq!(res.is_some(), label.starts_with("2F-pos"));
        insta::assert_yaml_snapshot!(res);

        Ok(())
    }

    /// Test internal workings of `handle_case_2h`.
    #[tracing_test::traced_test]
    #[rstest::rstest]
    // Note: same cases as in `evaluate` above -- keep in sync!
    // Case 2H: two or more HI predictors suggest AT LEAST ONE gene in the interval is HI
    #[case("1", 12_040_238, 12_073_572, "MFN2", "2H-neg")] // MFN2 not DECIPHER HI
    #[case("6", 111_620_236, 111_804_914, "REV3L", "2H-pos-1")] // REV3L, contained
    #[case("6", 111_620_236, 111_620_236, "REV3L", "2H-pos-2")] // REV3L, 1bp left
    #[case("6", 111_804_914, 111_804_914, "REV3L", "2H-pos-3")] // REV3L, 1bp right
    #[case("6", 111_620_235, 111_620_235, "REV3L", "2H-neg-1")] // REV3L, left of
    #[case("6", 111_804_915, 111_804_915, "REV3L", "2H-neg-2")] // REV3L, right of
    fn handle_case_2h(
        #[case] chrom: &str,
        #[case] start: u64,
        #[case] stop: u64,
        #[case] hgnc_id: &str,
        #[case] label: &str,
        global_evaluator_37: &super::super::super::Evaluator,
    ) -> Result<(), anyhow::Error> {
        mehari::common::set_snapshot_suffix!("{}-{}", hgnc_id, label);

        let strucvar = ds::StructuralVariant {
            chrom: chrom.into(),
            start: start as u32,
            stop: stop as u32,
            svtype: ds::SvType::Del,
            ambiguous_range: None,
        };
        let sv_interval = strucvar.into();

        let evaluator = super::Evaluator::with_parent(&global_evaluator_37);

        let res = evaluator.handle_case_2h(&sv_interval)?;
        insta::assert_yaml_snapshot!(res);

        Ok(())
    }

    /// Test internal workings of `get_gene_hi_predictions()`.`
    #[tracing_test::traced_test]
    #[rstest::rstest]
    fn get_gene_hi_predictions(
        global_evaluator_37: &super::super::super::Evaluator,
    ) -> Result<(), anyhow::Error> {
        let evaluator = super::Evaluator::with_parent(&global_evaluator_37);

        // MFN2 is not in DECIPHER HI.
        let res_mfn2 = evaluator.get_gene_hi_predictions("HGNC:16877");
        assert!(res_mfn2.is_none());

        // REV3L is in DECIPHER HI.
        let res_rev3l = evaluator.get_gene_hi_predictions("HGNC:9968");
        insta::assert_yaml_snapshot!(res_rev3l);

        // RERE is in DECIPHER HI.
        let res_rere = evaluator.get_gene_hi_predictions("HGNC:9965");
        insta::assert_yaml_snapshot!(res_rere);

        Ok(())
    }
}
