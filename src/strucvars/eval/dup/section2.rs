//! Implementation of evaluation of copy number gain section 2.

use bio::bio_types::genome::AbstractInterval;
use bio::bio_types::genome::Interval;

use crate::strucvars::eval::dup::result::G2C;
use crate::strucvars::{
    data::{
        clingen_dosage::{Gene, Region, Score},
        intervals::contains,
    },
    ds::StructuralVariant,
    eval::dup::result::{G2, G2B},
};

use super::result::{Section, G2A};

/// Evaluation of deletions, gain of copy number.
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

    /// Perform the evaluation of copy number gain Section 2 and all subsection.
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
    pub fn evaluate(&self, strucvar: &StructuralVariant) -> Result<Vec<Section>, anyhow::Error> {
        tracing::debug!("evaluating section 2 for {:?}", strucvar);
        let sv_interval: Interval = strucvar.clone().into();

        // Obtain any overlapping ClinGen region and gene dosage info.
        let (clingen_genes, clingen_regions) = self.parent.clingen_overlaps(&sv_interval);

        // Check each overlapping ClinGen region/gene.  If `sv_interval` completely contains
        // a region/gene and the region/gene has a "sufficient evidence" score for TS then we
        // have case 2A.
        if let Some(result) = Self::handle_case_2a(&clingen_genes, &clingen_regions, &sv_interval) {
            tracing::debug!("case 2A fired: {:?}", &result);
            return Ok(result);
        }
        // Otherwise, we need to test case 2B.
        let has_ts_genes = clingen_genes
            .iter()
            .any(|gene| gene.triplosensitivity_score == Score::SufficientEvidence);
        let has_ts_region = clingen_regions
            .iter()
            .any(|region| region.triplosensitivity_score == Score::SufficientEvidence);
        let mut result = if has_ts_genes || has_ts_region {
            let result = vec![Section::G2(G2::G2B(G2B::default()))];
            tracing::debug!("case 2B fired: {:?}", &result);
            result
        } else {
            tracing::debug!("negative for case 2B");
            Vec::new()
        };

        // Handle case 2C and return if it fired.
        let benign_regions = clingen_regions
            .iter()
            .filter(|region| region.triplosensitivity_score == Score::DosageSensitivityUnlikely)
            .collect::<Vec<_>>();
        if let Some(result_2c) = self.handle_case_2c(&sv_interval, &benign_regions)? {
            tracing::debug!("case 2C fired: {:?}", &result_2c);
            result.push(result_2c);
            return Ok(result);
        }

        Ok(result)
    }

    /// Handle case 2A.
    fn handle_case_2a(
        clingen_genes: &[Gene],
        clingen_regions: &[Region],
        sv_interval: &bio::bio_types::genome::Interval,
    ) -> Option<Vec<Section>> {
        let ts_genes = clingen_genes
            .iter()
            .filter(|clingen_gene| {
                let clingen_interval: Interval = (*clingen_gene)
                    .clone()
                    .try_into()
                    .expect("no interval for gene");
                contains(sv_interval, &clingen_interval)
                    && clingen_gene.triplosensitivity_score == Score::SufficientEvidence
            })
            .cloned()
            .collect::<Vec<_>>();
        let ts_regions = clingen_regions
            .iter()
            .cloned()
            .filter(|clingen_region| {
                let clingen_interval: Interval = clingen_region
                    .clone()
                    .try_into()
                    .expect("no interval for region");
                contains(sv_interval, &clingen_interval)
                    && clingen_region.triplosensitivity_score == Score::SufficientEvidence
            })
            .collect::<Vec<_>>();
        if !ts_regions.is_empty() || !ts_genes.is_empty() {
            Some(vec![Section::G2(G2::G2A(G2A {
                suggested_score: 1.0,
                ts_genes,
                ts_regions,
            }))])
        } else {
            None
        }
    }

    // Handle case 2C (identical in gene content to the established benign copy-number gain).
    fn handle_case_2c(
        &self,
        sv_interval: &Interval,
        benign_regions: &[&Region],
    ) -> Result<Option<Section>, anyhow::Error> {
        // Short-circuit if there are no TS regions.
        if benign_regions.is_empty() {
            return Ok(None);
        }

        // Get sorted HGNC IDs of overlapping genes.
        let sv_hgncs = self.parent.overlapping_gene_hcnc_ids(
            sv_interval.contig(),
            sv_interval.range().start as u32,
            sv_interval.range().end as u32,
        )?;

        // Now, check each overlapping TS region for identical gene content.
        let mut res_benign_regions = Vec::new();
        for benign_region in benign_regions {
            let ts_interval: Interval = (*benign_region).clone().try_into().map_err(|e| {
                anyhow::anyhow!("failed to convert genomic location to interval: {}", e)
            })?;
            let region_hgncs = self.parent.overlapping_gene_hcnc_ids(
                ts_interval.contig(),
                ts_interval.range().start as u32,
                ts_interval.range().end as u32,
            )?;
            if sv_hgncs == region_hgncs {
                res_benign_regions.push((*benign_region).clone());
            }
        }

        if res_benign_regions.is_empty() {
            Ok(None)
        } else {
            Ok(Some(Section::G2(G2::G2C(G2C {
                suggested_score: -1.0,
                benign_regions: res_benign_regions,
            }))))
        }
    }
}

#[cfg(test)]
pub mod test {
    use hgvs::data::interface::Provider;

    use crate::strucvars::data::clingen_dosage::Score;
    use crate::strucvars::ds;

    use super::super::super::test::global_evaluator_37;

    /// Test internal working of `handle_case_2a` (complete overlap).
    #[tracing_test::traced_test]
    #[rstest::rstest]
    // Note: same cases as in `evaluate` above -- keep in sync!
    #[case("X", 123_034_319, 123_236_519, "ISCA-46743", "2A-pos")]
    #[case("X", 103_031_434, 103_047_548, "LMNB1", "2A-pos")]
    fn handle_case_2a(
        #[case] chrom: &str,
        #[case] start: u64,
        #[case] stop: u64,
        #[case] hgnc_id: &str,
        #[case] label: &str,
        global_evaluator_37: super::super::super::Evaluator,
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

    /// Test internal working of `handle_case_2c` (identical in gene content with benign).
    #[tracing_test::traced_test]
    #[rstest::rstest]
    // Note: same cases as in `evaluate` above -- keep in sync!
    #[case("2", 110_862_108, 110_983_703, "ISCA-37405", "2C-pos-1")] // exact match
    #[case("2", 110_802_355, 110_984_922, "ISCA-37405", "2C-pos-2")] // full: MALL, MTLN
    #[case("2", 110_862_108, 111_246_742, "ISCA-37405", "2C-neg-1")] // extra: LIMS4
    #[case("2", 110_879_569, 110_983_703, "ISCA-37405", "2C-neg-2")] // excludes: MALL
    fn handle_case_2c(
        #[case] chrom: &str,
        #[case] start: u64,
        #[case] stop: u64,
        #[case] hgnc_id: &str,
        #[case] label: &str,
        global_evaluator_37: super::super::super::Evaluator,
    ) -> Result<(), anyhow::Error> {
        mehari::common::set_snapshot_suffix!("{}-{}", hgnc_id, label);
        let evaluator = super::Evaluator::with_parent(&global_evaluator_37);

        let strucvar = ds::StructuralVariant {
            chrom: chrom.into(),
            start: start as u32,
            stop: stop as u32,
            svtype: ds::SvType::Del,
            ambiguous_range: None,
        };
        let sv_interval = strucvar.into();

        let (_, clingen_regions) = global_evaluator_37.clingen_overlaps(&sv_interval);
        let benign_regions = clingen_regions
            .iter()
            .filter(|region| region.triplosensitivity_score == Score::DosageSensitivityUnlikely)
            .collect::<Vec<_>>();

        let res = evaluator.handle_case_2c(&sv_interval, &benign_regions)?;
        insta::assert_yaml_snapshot!(res);

        Ok(())
    }
}
