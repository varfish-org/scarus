//! Implementation of evaluation of copy number loss section 4.
//!
//! Note that we can only determine 4O reliably in an automated fashion.  We support
//! reporting overlapping variants from ClinVar but leave them as "dangle" because a
//! human must evaluate the phenotype.

use super::result::{L4Patho, Section, L4, L4N};
use crate::strucvars::{
    ds::StructuralVariant,
    eval::{
        del::result::L4O,
        result::{ClinvarSvOverlap, GnomadSvOverlap},
        CarrierFrequency, IntoInterval, Overlaps, SvType,
    },
};
use annonars::{
    clinvar_sv::cli::query::EXAC_CNV_CASES,
    gnomad_sv::cli::query::Record as GnomadSvRecord,
    pbs::annonars::clinvar::v1::{
        minimal::{ClinicalSignificance, ReviewStatus, VariantType},
        sv::Record as ClinvarSvRecord,
    },
};
use bio::bio_types::genome::{AbstractInterval, Interval};

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

    /// Perform the evaluation of copy number loss Section 4.
    ///
    /// Note that the results depend on the phenotype of the patient when considering ClinVar
    /// SVs.  We thus report pathogenic and benign variants in ClinVar and gnomAD separately
    /// and a human must resolve them.
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
        tracing::debug!("evaluating section 4 for {:?}", strucvar);
        let sv_interval: Interval = strucvar.clone().into();
        let mut result = Vec::new();

        // Retrieve all overlapping ClinVar records.
        let clinvar_records = self.overlapping_clinvar_records(&sv_interval)?;

        if let Some(result_clinvar_patho) =
            self.handle_clinvar_pathogenic(&sv_interval, &clinvar_records)?
        {
            result.push(result_clinvar_patho)
        }
        if let Some(result_clinvar_benign) =
            self.handle_clinvar_benign(&sv_interval, &clinvar_records)?
        {
            result.push(result_clinvar_benign)
        }
        if let Some(result_gnomad_sv) = self.handle_gnomad_sv(&sv_interval)? {
            result.push(result_gnomad_sv)
        }

        tracing::warn!("Section 4 evaluation not implemented yet");
        Ok(result)
    }

    /// Obtain overlapping ClinVar records.
    fn overlapping_clinvar_records(
        &self,
        sv_interval: &Interval,
    ) -> Result<Vec<ClinvarSvRecord>, anyhow::Error> {
        tracing::debug!("querying ClinVar for {:?}", sv_interval);
        let res = self
            .parent
            .clinvar_sv_overlaps(
                sv_interval.contig(),
                sv_interval.range().start as u32,
                sv_interval.range().end as u32,
            )
            .map_err(|e| anyhow::anyhow!("ClinVar query failed: {}", e))?;
        tracing::debug!("found {} overlapping ClinVar records", res.len());
        tracing::trace!("records are {:#?}", &res);
        Ok(res)
    }

    /// Evaluate ClinVar pathogenic variants.
    ///
    /// This could be one of 4A, 4B, 4C, 4D, 4E, 4L, 4M.
    fn handle_clinvar_pathogenic(
        &self,
        sv_interval: &Interval,
        clinvar_records: &[ClinvarSvRecord],
    ) -> Result<Option<Section>, anyhow::Error> {
        let pathogenic = [
            ClinicalSignificance::Pathogenic as i32,
            ClinicalSignificance::LikelyPathogenic as i32,
        ];
        // Filter intervals to pathogenic ones without conflicts and sufficient overlap.
        let mut overlaps = clinvar_records
            .iter()
            .flat_map(|record| {
                let rc = record.reference_assertions.first().expect("no RCV");
                let record_interval = record.into_interval();
                let overlap = sv_interval.reciprocal_overlap(&record_interval) as f32;
                tracing::trace!("considering records {:?} with overlap {}", record, overlap);
                tracing::trace!("rc.clinical_significance = {} / {:?}", rc.clinical_significance, &pathogenic);
                let is_pathogenic = pathogenic.contains(&rc.clinical_significance);
                let is_review_status = rc.review_status <= ReviewStatus::CriteriaProvidedSingleSubmitter as i32;
                let is_variant_type = record.variant_type == VariantType::Deletion as i32;
                let is_overlap = overlap >= self.parent.config.clinvar_sv_min_overlap_pathogenic;
                let is_interesting = is_pathogenic && is_variant_type && is_overlap;
                tracing::trace!("is_pathogenic: {}, is_review_status: {}, is_variant_type: {}, is_overlap: {}, is_interesting: {}", is_pathogenic, is_review_status, is_variant_type, is_overlap, is_interesting);
                if is_interesting {
                    Some(ClinvarSvOverlap {
                        overlap,
                        record: record.clone(),
                    })
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();

        overlaps.sort_by(|a, b| {
            b.overlap
                .partial_cmp(&a.overlap)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        if overlaps.is_empty() {
            Ok(None)
        } else {
            Ok(Some(Section::L4(L4::L4Patho(L4Patho { overlaps }))))
        }
    }

    /// Evaluate ClinVar benign variants.
    ///
    /// This is roughly 4N.
    fn handle_clinvar_benign(
        &self,
        sv_interval: &Interval,
        clinvar_records: &[ClinvarSvRecord],
    ) -> Result<Option<Section>, anyhow::Error> {
        let benign = [
            ClinicalSignificance::Benign as i32,
            ClinicalSignificance::LikelyBenign as i32,
        ];
        // Filter intervals to pathogenic ones without conflicts and sufficient overlap.
        let mut overlaps = clinvar_records
            .iter()
            .flat_map(|record| {
                let rc = record.reference_assertions.first().expect("no RCV");
                let record_interval = record.into_interval();
                let overlap = sv_interval.reciprocal_overlap(&record_interval) as f32;
                tracing::trace!("considering records {:?} with overlap {}", record, overlap);
                let interesting = benign.contains(&rc.clinical_significance)
                    && rc.review_status <= ReviewStatus::CriteriaProvidedSingleSubmitter as i32
                    && record.variant_type == VariantType::Deletion as i32
                    && overlap >= self.parent.config.clinvar_sv_min_overlap_pathogenic;
                if interesting {
                    Some(ClinvarSvOverlap {
                        overlap,
                        record: record.clone(),
                    })
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();

        overlaps.sort_by(|a, b| {
            b.overlap
                .partial_cmp(&a.overlap)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        if overlaps.is_empty() {
            Ok(None)
        } else {
            Ok(Some(Section::L4(L4::L4N(L4N { overlaps }))))
        }
    }

    /// Handle overlap in gnomAD-SV.
    ///
    /// This corresponds to case 4O.
    fn handle_gnomad_sv(&self, sv_interval: &Interval) -> Result<Option<Section>, anyhow::Error> {
        // Query gnomAD-SV for all records.
        tracing::debug!("querying gnomAD-SV for {:?}", sv_interval);
        let raw_records = self
            .parent
            .gnomad_sv_overlaps(
                sv_interval.contig(),
                sv_interval.range().start as u32,
                sv_interval.range().end as u32,
            )
            .map_err(|e| anyhow::anyhow!("gnomAD-SV query failed: {}", e))?;
        tracing::debug!("found {} overlapping gnomAD-SV records", raw_records.len());
        tracing::trace!("records are {:#?}", &raw_records);

        // Filter to records with the correct type and appropriate overlap.
        let mut overlaps: Vec<GnomadSvOverlap> = raw_records
            .into_iter()
            .flat_map(|record| {
                let record_interval = record.into_interval();
                let overlap = sv_interval.reciprocal_overlap(&record_interval) as f32;
                let sv_type = SvType::try_from(&record);
                let interesting = matches!(sv_type, Ok(SvType::Del | SvType::Cnv))
                    && overlap >= self.parent.config.gnomad_sv_min_overlap;
                if interesting {
                    Some(GnomadSvOverlap { overlap, record })
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();
        overlaps.sort_by(|a, b| {
            b.overlap
                .partial_cmp(&a.overlap)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        tracing::trace!("filtered overlaps are {:#?}", &overlaps);

        // We have to manually count ExAC CNV records as they are per-sample and not aggregated.
        let exac_cnv_count = overlaps
            .iter()
            .filter(|GnomadSvOverlap { record, .. }| matches!(record, &GnomadSvRecord::ExacCnv(_)))
            .count();
        let exac_freq = exac_cnv_count as f32 / EXAC_CNV_CASES as f32;
        // For the rest, we can obtain the frequency from the records.
        let other_freqs = overlaps
            .iter()
            .flat_map(|GnomadSvOverlap { record, .. }| {
                if matches!(record, &GnomadSvRecord::ExacCnv(_)) {
                    None
                } else {
                    Some(record.carrier_freq().expect("no carrier freq?"))
                }
            })
            .collect::<Vec<_>>();
        let max_freq = other_freqs
            .iter()
            .max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .copied()
            .unwrap_or(0.0);

        let suggested_score = if exac_freq >= self.parent.config.gnomad_sv_min_freq_common
            || max_freq >= self.parent.config.gnomad_sv_min_freq_common
        {
            -1.0
        } else {
            0.0
        };

        Ok(if overlaps.is_empty() {
            None
        } else {
            Some(Section::L4(L4::L4O(L4O {
                suggested_score,
                overlaps,
            })))
        })
    }
}

#[cfg(test)]
mod test {
    use crate::strucvars::ds::{StructuralVariant, SvType};

    use super::super::super::test::global_evaluator_37;

    #[tracing_test::traced_test]
    #[rstest::rstest]
    // ClinVar pathogenic variant
    #[case("17", 43_005_866, 43_377_096, "VCV000059587")]
    // ClinVar benign variant
    #[case("X", 155_231_256, 155_251_871, "VCV000161054")]
    // gnomAD SV variants
    // ExAC CNV
    #[case("10", 51_620_321, 51_748_701, "chr10-51620321-51748701-NFE")]
    // gnomAD-SV v2
    #[case("1", 21000, 26000, "gnomAD-SV_v2.1_DEL_1_")]
    fn test_evaluate(
        #[case] chrom: &str,
        #[case] start: u32,
        #[case] stop: u32,
        #[case] label: &str,
        global_evaluator_37: &super::super::super::Evaluator,
    ) -> Result<(), anyhow::Error> {
        mehari::common::set_snapshot_suffix!("{}", label);

        let evaluator = super::Evaluator::with_parent(global_evaluator_37);

        let strucvar = StructuralVariant {
            chrom: chrom.into(),
            start: start,
            stop: stop,
            svtype: SvType::Del,
            ambiguous_range: None,
        };

        insta::assert_yaml_snapshot!(evaluator.evaluate(&strucvar)?);

        Ok(())
    }

    #[tracing_test::traced_test]
    #[rstest::rstest]
    // same as in `evaluate` test above, keep in sync
    #[case("17", 43_005_866, 43_377_096, "VCV000059587")]
    fn handle_clinvar_pathogenic(
        #[case] chrom: &str,
        #[case] start: u32,
        #[case] stop: u32,
        #[case] label: &str,
        global_evaluator_37: &super::super::super::Evaluator,
    ) -> Result<(), anyhow::Error> {
        mehari::common::set_snapshot_suffix!("{}", label);

        let evaluator = super::Evaluator::with_parent(global_evaluator_37);
        let strucvar = StructuralVariant {
            chrom: chrom.into(),
            start: start,
            stop: stop,
            svtype: SvType::Del,
            ambiguous_range: None,
        };
        let sv_interval = strucvar.into();

        let clinvar_records = evaluator.overlapping_clinvar_records(&sv_interval)?;

        insta::assert_yaml_snapshot!(
            evaluator.handle_clinvar_pathogenic(&sv_interval, &clinvar_records)?
        );

        Ok(())
    }

    #[tracing_test::traced_test]
    #[rstest::rstest]
    // same as in `evaluate` test above, keep in sync
    #[case("X", 155_231_256, 155_251_871, "VCV000161054")]
    fn handle_clinvar_benign(
        #[case] chrom: &str,
        #[case] start: u32,
        #[case] stop: u32,
        #[case] label: &str,
        global_evaluator_37: &super::super::super::Evaluator,
    ) -> Result<(), anyhow::Error> {
        mehari::common::set_snapshot_suffix!("{}", label);

        let evaluator = super::Evaluator::with_parent(global_evaluator_37);
        let strucvar = StructuralVariant {
            chrom: chrom.into(),
            start: start,
            stop: stop,
            svtype: SvType::Del,
            ambiguous_range: None,
        };
        let sv_interval = strucvar.into();

        let clinvar_records = evaluator.overlapping_clinvar_records(&sv_interval)?;

        insta::assert_yaml_snapshot!(
            evaluator.handle_clinvar_benign(&sv_interval, &clinvar_records)?
        );

        Ok(())
    }

    #[tracing_test::traced_test]
    #[rstest::rstest]
    // same as in `evaluate` test above, keep in sync
    // ExAC CNV
    #[case("10", 51_620_321, 51_748_701, "chr10-51620321-51748701-NFE")]
    // gnomAD-SV v2
    #[case("1", 21000, 26000, "gnomAD-SV_v2.1_DEL_1_")]
    fn handle_gnomad_sv(
        #[case] chrom: &str,
        #[case] start: u32,
        #[case] stop: u32,
        #[case] label: &str,
        global_evaluator_37: &super::super::super::Evaluator,
    ) -> Result<(), anyhow::Error> {
        mehari::common::set_snapshot_suffix!("{}", label);

        let evaluator = super::Evaluator::with_parent(global_evaluator_37);
        let strucvar = StructuralVariant {
            chrom: chrom.into(),
            start: start,
            stop: stop,
            svtype: SvType::Del,
            ambiguous_range: None,
        };
        let sv_interval = strucvar.into();

        insta::assert_yaml_snapshot!(evaluator.handle_gnomad_sv(&sv_interval)?);

        Ok(())
    }
}
