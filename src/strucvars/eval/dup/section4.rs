//! Implementation of evaluation of copy number loss section 4.
//!
//! Note that we can only determine 4O reliably in an automated fashion.  We support
//! reporting overlapping variants from ClinVar but leave them as "dangle" because a
//! human must evaluate the phenotype.

use super::result::{G4Patho, Section, G4, G4N};
use crate::strucvars::{
    ds::StructuralVariant,
    eval::{
        dup::result::G4O,
        result::{ClinvarSvOverlap, GnomadSvOverlap},
        CarrierFrequency, IntoInterval, Overlaps, SvType,
    },
};
use annonars::{
    clinvar_sv::cli::query::EXAC_CNV_CASES,
    gnomad_sv::cli::query::Record as GnomadSvRecord,
    pbs::clinvar::{
        minimal::{ClinicalSignificance, ReviewStatus, VariantType},
        sv::Record as ClinvarSvRecord,
    },
};
use bio::bio_types::genome::{AbstractInterval, Interval};

/// Evaluation of duplications, gain of copy number.
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
                let is_variant_type = record.variant_type == VariantType::Duplication as i32;
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
            Ok(Some(Section::G4(G4::G4Patho(G4Patho { overlaps }))))
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
                    && record.variant_type == VariantType::Duplication as i32
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
            Ok(Some(Section::G4(G4::G4N(G4N { overlaps }))))
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
                let interesting = matches!(sv_type, Ok(SvType::Dup | SvType::Cnv))
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
            Some(Section::G4(G4::G4O(G4O {
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
    #[case("10", 55_581_618, 56_288_162, "VCV000806586")]
    // ClinVar benign variant
    #[case("10", 135_252_347, 135_378_802, "VCV000221841")]
    // gnomAD SV variants
    // ExAC CNV
    #[case("10", 92_994, 95_183, "chr10-92994-95183-NFE")]
    // gnomAD-SV v2
    #[case("1", 157_000, 166_000, "gnomAD-SV_v2.1_DUP_1_6")]
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
            start,
            stop,
            svtype: SvType::Dup,
            ambiguous_range: None,
        };

        insta::assert_yaml_snapshot!(evaluator.evaluate(&strucvar)?);

        Ok(())
    }

    #[tracing_test::traced_test]
    #[rstest::rstest]
    // same as in `evaluate` test above, keep in sync
    #[case("10", 55_581_618, 56_288_162, "VCV000806586")]
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
            start,
            stop,
            svtype: SvType::Dup,
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
    #[case("10", 135_252_347, 135_378_802, "VCV000221841")]
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
            start,
            stop,
            svtype: SvType::Dup,
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
    #[case("10", 92_994, 95_183, "chr10-92994-95183-NFE")]
    // gnomAD-SV v2
    #[case("1", 157_000, 166_000, "gnomAD-SV_v2.1_DUP_1_6")]
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
            start,
            stop,
            svtype: SvType::Dup,
            ambiguous_range: None,
        };
        let sv_interval = strucvar.into();

        insta::assert_yaml_snapshot!(evaluator.handle_gnomad_sv(&sv_interval)?);

        Ok(())
    }
}
