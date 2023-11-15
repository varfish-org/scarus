//! Implementation of evaluation of copy number gain section 2.

use bio::bio_types::genome::AbstractInterval;
use bio::bio_types::genome::Interval;
use hgvs::data::interface::Provider as _;

use crate::strucvars::eval::dup::result::G2L;
use crate::strucvars::eval::result::Pvs1Result;
use crate::strucvars::{
    data::{
        clingen_dosage::{Gene, Region, Score},
        hgnc::GeneIdInfo,
        intervals::contains,
        intervals::do_overlap,
    },
    ds::StructuralVariant,
    eval::dup::result::{G2, G2B, G2C, G2D, G2E, G2F, G2G, G2H, G2I, G2J},
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

        // Handle cases 2D..2G and return if 2D/2F fired.
        let result_2d_2g = self.handle_cases_2d_2g(&sv_interval, &benign_regions)?;
        for entry in &result_2d_2g {
            tracing::debug!("case 2D..2G fired: {:?}", &entry);
            match entry {
                Section::G2(G2::G2D(_)) | Section::G2(G2::G2F(_)) => {
                    tracing::debug!("case 2D/2F fired");
                    result.push(entry.clone());
                    return Ok(result);
                }
                Section::G2(G2::G2E(_)) | Section::G2(G2::G2G(_)) => {
                    tracing::debug!("case 2E/2G fired");
                    result.push(entry.clone());
                    // continue evaluating
                }
                _ => unreachable!(),
            }
        }

        // Handle cases 2H..2L.
        let mut result_2h_2l = self.handle_cases_2h_2l(&sv_interval, &clingen_genes)?;
        result.append(&mut result_2h_2l);
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

    /// Handle case 2C (identical in gene content to the established benign copy-number gain).
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

    /// Handle cases 2D..2G.
    ///
    /// Precondition: negative for 2C, not identical in gene content to established benign region.
    fn handle_cases_2d_2g(
        &self,
        sv_interval: &Interval,
        benign_regions: &[&Region],
    ) -> Result<Vec<Section>, anyhow::Error> {
        // Get overlapping genes and transcripts.
        let sv_genes = self.parent.overlapping_genes(
            sv_interval.contig(),
            sv_interval.range().start as u32,
            sv_interval.range().end as u32,
        )?;
        let contig_ac = self
            .parent
            .chrom_to_ac
            .get(sv_interval.contig())
            .expect("unknown contig");

        // Obtain CDSs of overlapping coding transcripts.
        let coding_tx_cds = sv_genes
            .iter()
            .flat_map(|overlap| {
                overlap
                    .tx_ids
                    .iter()
                    .flat_map(|tx_id| -> Option<(GeneIdInfo, Interval)> {
                        self.parent.provider.get_tx(tx_id).and_then(
                            |tx| -> Option<(GeneIdInfo, Interval)> {
                                for genome_alignment in &tx.genome_alignments {
                                    if let (Some(cds_start), Some(cds_end)) =
                                        (genome_alignment.cds_start, genome_alignment.cds_end)
                                    {
                                        if contig_ac == &genome_alignment.contig {
                                            return Some((
                                                overlap.gene.clone(),
                                                Interval::new(
                                                    sv_interval.contig().to_string(),
                                                    (cds_start as u64)..(cds_end as u64),
                                                ),
                                            ));
                                        }
                                    }
                                }
                                None
                            },
                        )
                    })
            })
            .collect::<Vec<_>>();

        // Cases 2D/2E: SV is smaller than the benign interval
        //
        // Collect potentially interrupted gene HGNC IDs.
        let mut containing_region = None;
        let mut result = Vec::new();
        for &benign_region in benign_regions {
            let benign_interval: Interval = benign_region.clone().try_into().map_err(|e| {
                anyhow::anyhow!("failed to convert genomic location to interval: {}", e)
            })?;

            if contains(&benign_interval, sv_interval) {
                containing_region = Some(benign_region.clone());
                let mut interrupted_genes = Vec::new();
                for (hgnc_id, cds) in &coding_tx_cds {
                    // Any breakpoint of sv_interval within CDS is the same as
                    // having an overlap but none containg the CDS.
                    if do_overlap(sv_interval, cds) && !contains(sv_interval, cds) {
                        interrupted_genes.push(hgnc_id.clone());
                    }
                }
                interrupted_genes.sort_by(|a, b| a.hgnc_id.cmp(&b.hgnc_id));
                interrupted_genes.dedup_by(|a, b| a.hgnc_id == b.hgnc_id);

                if !interrupted_genes.is_empty() {
                    tracing::debug!("case 2E fired: {:?}", &interrupted_genes);
                    // Case 2E, potentially interrupts coding region.
                    result.push(Section::G2(G2::G2E(G2E {
                        suggested_score: 0.0,
                        benign_region: benign_region.clone(),
                        genes: interrupted_genes,
                    })));
                }
            }
        }
        if containing_region.is_some() && result.is_empty() {
            tracing::debug!("case 2D fired: {:?}", &containing_region);
            return Ok(vec![Section::G2(G2::G2D(G2D {
                suggested_score: -1.0,
                benign_region: containing_region.expect("checked above"),
            }))]);
        }

        // If we reach here, then the DUP is larger than the region.
        let mut contained_region = None; // for 2F
        let mut found_additional = false;
        for &benign_region in benign_regions {
            let benign_interval: Interval = benign_region.clone().try_into().map_err(|e| {
                anyhow::anyhow!("failed to convert genomic location to interval: {}", e)
            })?;

            // Get all overlapping genes in benign region.
            let benign_genes = self
                .parent
                .overlapping_gene_hcnc_ids(
                    benign_interval.contig(),
                    benign_interval.range().start as u32,
                    benign_interval.range().end as u32,
                )?
                .into_iter()
                .collect::<rustc_hash::FxHashSet<_>>();
            // Check whether all protein-coding genes of the SV are also in the benign region.
            let mut additional_genes: Vec<GeneIdInfo> = coding_tx_cds
                .iter()
                .flat_map(|(gene_id_info, _)| {
                    if benign_genes.contains(&gene_id_info.hgnc_id) {
                        None
                    } else {
                        Some(gene_id_info.clone())
                    }
                })
                .collect::<Vec<_>>();
            additional_genes.sort_by(|a, b| a.hgnc_id.cmp(&b.hgnc_id));
            additional_genes.dedup_by(|a, b| a.hgnc_id == b.hgnc_id);

            if additional_genes.is_empty() {
                // For 2F, SV must be larger than benign interval (not the same).
                tracing::trace!(
                    "sv_interval = {:?}, benign_interval = {:?}",
                    sv_interval,
                    &benign_interval
                );
                if contains(sv_interval, &benign_interval)
                    && !contains(&benign_interval, sv_interval)
                {
                    contained_region = Some(benign_region.clone());
                }
            } else {
                tracing::debug!("case 2G fired: {:?}", &additional_genes);
                found_additional = true;
                result.push(Section::G2(G2::G2G(G2G {
                    suggested_score: 0.0,
                    benign_region: benign_region.clone(),
                    genes: additional_genes,
                })));
            }
        }
        tracing::trace!(
            "found_additional: {}, contained_region: {:?}",
            found_additional,
            &contained_region
        );
        if !found_additional {
            if let Some(contained_region) = contained_region {
                tracing::debug!("case 2F fired");
                result.push(Section::G2(G2::G2F(G2F {
                    suggested_score: -1.0,
                    benign_region: contained_region,
                })));
            }
        }

        Ok(result)
    }

    /// Handle cases 2H..2L.
    ///
    /// Note that we cannot reiliably decide for 2K, so we will only return 2J.
    fn handle_cases_2h_2l(
        &self,
        sv_interval: &Interval,
        clingen_genes: &[Gene],
    ) -> Result<Vec<Section>, anyhow::Error> {
        let mut result = Vec::new();

        {
            // For case 2H: get HGNC IDs of overlapping HI genes
            let hi_genes = clingen_genes
                .iter()
                .filter(|gene| {
                    let gene_interval: Interval =
                        (*gene).clone().try_into().expect("conversion failed");
                    gene.haploinsufficiency_score == Score::SufficientEvidence
                        && contains(sv_interval, &gene_interval)
                })
                .flat_map(|gene| {
                    self.parent
                        .gene_id_data
                        .by_ncbi_gene_id(&gene.ncbi_gene_id)
                        .cloned()
                })
                .collect::<Vec<_>>();
            if !hi_genes.is_empty() {
                // Case 2H: HI gene(s) are completely contained in the SV.
                tracing::debug!("case 2H fired: {:?}", &hi_genes);
                result.push(Section::G2(G2::G2H(G2H {
                    suggested_score: 0.0,
                    hi_genes,
                })));
            }
        }

        {
            // For case 2I: both breakpoints in the same HI gene.
            let hi_genes = clingen_genes
                .iter()
                .filter(|gene| {
                    let gene_interval: Interval =
                        (*gene).clone().try_into().expect("conversion failed");
                    // no reciprocal containment => true containment
                    gene.haploinsufficiency_score == Score::SufficientEvidence
                        && contains(&gene_interval, sv_interval)
                        && !contains(sv_interval, &gene_interval)
                })
                .flat_map(|gene| {
                    self.parent
                        .gene_id_data
                        .by_ncbi_gene_id(&gene.ncbi_gene_id)
                        .cloned()
                })
                .collect::<Vec<_>>();
            if !hi_genes.is_empty() {
                // Case 2I: HI gene(s) are completely contained in the SV.
                tracing::debug!("case 2I fired: {:?}", &hi_genes);
                // TODO: proper implementat with PVS1
                result.push(Section::G2(G2::G2I(G2I {
                    suggested_score: 0.9,
                    pvs1_result: Pvs1Result::Pvs1,
                    hi_genes,
                })));
                return Ok(result);
            }
        }

        {
            // For case 2K: one breakpoint in HI gene.
            //
            // As mentioned above, we won't do 2J as we cannot incorporate the phenotype.
            let hi_genes = clingen_genes
                .iter()
                .filter(|gene| {
                    let gene_interval: Interval =
                        (*gene).clone().try_into().expect("conversion failed");
                    // We look for true overlap here, no equality
                    gene.haploinsufficiency_score == Score::SufficientEvidence
                        && do_overlap(sv_interval, &gene_interval)
                        && !contains(sv_interval, &gene_interval)
                        && !contains(&gene_interval, sv_interval)
                })
                .flat_map(|gene| {
                    self.parent
                        .gene_id_data
                        .by_ncbi_gene_id(&gene.ncbi_gene_id)
                        .cloned()
                })
                .collect::<Vec<_>>();
            if !hi_genes.is_empty() {
                // Case 2J: one breakpoint is in HI gene
                tracing::debug!("case 2J fired: {:?}", &hi_genes);
                result.push(Section::G2(G2::G2J(G2J {
                    suggested_score: 0.45,
                    hi_genes,
                })));
                return Ok(result);
            }
        }

        {
            // For case 2L: one breakpoint in any gene.

            // Map chromosome name (e.g., chr1) to chromosome accession in this assembly.
            let chrom_acc = self
                .parent
                .chrom_to_ac
                .get(sv_interval.contig())
                .ok_or_else(|| {
                    anyhow::anyhow!("could not resolve contig name `{}`", sv_interval.contig())
                })?;

            // Obtain the overlapping transcripts for the given region.
            let txs = self
                .parent
                .provider
                .get_tx_for_region(
                    chrom_acc,
                    "splign",
                    sv_interval.range().start as i32,
                    sv_interval.range().end as i32,
                )
                .map_err(|e| {
                    anyhow::anyhow!("problem query transcript database with range: {}", e)
                })?
                .iter()
                .map(|tx| {
                    self.parent
                        .provider
                        .get_tx_exons(&tx.tx_ac, chrom_acc, "splign")
                })
                .collect::<Result<Vec<_>, _>>()?;

            tracing::trace!("txs = {:?}", &txs);

            // Obtain truly overlapping genes.
            let mut hgnc_ids = Vec::new();
            for tx_exons in &txs {
                let tx_start = tx_exons
                    .iter()
                    .map(|exon| exon.alt_start_i - 1)
                    .min()
                    .expect("no exons?") as u64;
                let tx_end = tx_exons
                    .iter()
                    .map(|exon| exon.alt_end_i)
                    .max()
                    .expect("no exons?") as u64;
                let tx_interval = Interval::new(sv_interval.contig().to_string(), tx_start..tx_end);
                if !contains(sv_interval, &tx_interval) && !contains(&tx_interval, sv_interval) {
                    hgnc_ids.push(
                        self.parent
                            .provider
                            .get_tx_identity_info(&tx_exons.first().expect("no exons?").tx_ac)?
                            .hgnc,
                    );
                }
            }
            hgnc_ids.sort();
            hgnc_ids.dedup();

            if !hgnc_ids.is_empty() {
                let genes = hgnc_ids
                    .iter()
                    .flat_map(|hgnc_id| self.parent.gene_id_data.by_hgnc_id(hgnc_id).cloned())
                    .collect::<Vec<_>>();

                // Case 2L: one breakpoint is in any gene
                tracing::debug!("case 2L fired: {:?}", &hgnc_ids);
                result.push(Section::G2(G2::G2L(G2L {
                    suggested_score: 0.45,
                    genes,
                })));
            }
        }

        Ok(result)
    }
}

#[cfg(test)]
pub mod test {
    use crate::strucvars::data::clingen_dosage::Score;
    use crate::strucvars::ds;

    use super::super::super::test::global_evaluator_37;

    /// Test the `evaluate` function with the different cases of Section 2.
    #[tracing_test::traced_test]
    #[rstest::rstest]
    // Case 2A
    #[case("X", 123_034_319, 123_236_519, "ISCA-46743", "2A-pos")]
    #[case("X", 103_031_434, 103_047_548, "LMNB1", "2A-pos")]
    // Case 2C
    #[case("2", 110_862_108, 110_983_703, "ISCA-37405", "2C-pos-1")] // exact match
    #[case("2", 110_802_355, 110_984_922, "ISCA-37405", "2C-pos-2")] // full: MALL, MTLN
    #[case("2", 110_862_108, 111_246_742, "ISCA-37405", "2C-neg-1")] // extra: LIMS4
    #[case("2", 110_879_569, 110_983_703, "ISCA-37405", "2C-neg-2")]
    // excludes: MALL
    // Case 2D: smaller than established benign region, no coding genes interrupted.
    #[case("2", 110_876_929, 110_965_509, "ISCA-37405", "2D-pos-1")]
    // smaller, no interrupt
    // Case 2E: smaller than established benign region, potential protein coding interrupt.
    #[case("2", 110_862_109, 110_983_702, "ISCA-37405", "2E-pos-1")] // region (-1bp) interrupts
    #[case("2", 110_878_786, 110_969_992, "ISCA-37405", "2E-pos-2")] // interrupt MTLN
    #[case("2", 110_882_589, 110_961_118, "ISCA-37405", "2E-pos-3")]
    // interrupt NPHP1
    // Case 2F: larger than established benign region, no additional genetic material.
    #[case("2", 110_862_107, 110_983_704, "ISCA-37405", "2F-pos-1")] // +1bp
    // Case 2G: larger than established benign region, additional protein-coding gene.
    #[case("2", 110_833_712, 111_236_476, "ISCA-37405", "2G-pos-1")] // adds LIMS4
    // Case 2H: HI gene fully contained within observed copy number gain.
    #[case("2", 189839099, 189877472, "COL3A1", "2H-pos-1")]
    // Case 2I: Both breakpoints in the same HI gene.
    #[case("2", 189839100, 189877471, "COL3A1", "2I-pos-1")]
    // Case 2J: One breakpoint in an HI gene.
    #[case("2", 189877072, 189877500, "COL3A1", "2J-pos-1")] // left bp
    #[case("2", 189839000, 189839599, "COL3A1", "2J-pos-2")]
    // right bp
    // No Case 2K: we cannot handle phenotypes well.
    // Case 2L: One breakpoint in any gene.
    #[case("6", 111_620_000, 111_621_000, "REV3L", "2L-pos-1")] // REV3L left bp
    #[case("6", 111_800_000, 111_805_000, "REV3L", "2L-pos-1")] // REV3L right bp
    fn evaluate(
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

        let res = evaluator.evaluate(&strucvar)?;
        insta::assert_yaml_snapshot!(res);

        Ok(())
    }

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

    /// Test internal working of `handle_cases_2d_2g`.
    #[tracing_test::traced_test]
    #[rstest::rstest]
    // Note: same cases as in `evaluate` above -- keep in sync!
    //
    // Case 2D: smaller than established benign region, no coding genes interrupted.
    #[case("2", 110_876_929, 110_965_509, "ISCA-37405", "2D-pos-1")]
    // smaller, no interrupt
    // Case 2E: smaller than established benign region, potential protein coding interrupt.
    #[case("2", 110_862_109, 110_983_702, "ISCA-37405", "2E-pos-1")] // region (-1bp) interrupts
    #[case("2", 110_878_786, 110_969_992, "ISCA-37405", "2E-pos-2")] // interrupt MTLN
    #[case("2", 110_882_589, 110_961_118, "ISCA-37405", "2E-pos-3")]
    // interrupt NPHP1
    // Case 2F: larger than established benign region, no additional genetic material.
    #[case("2", 110_862_107, 110_983_704, "ISCA-37405", "2F-pos-1")] // +1bp
    // Case 2G: larger than established benign region, additional protein-coding gene.
    #[case("2", 110_833_712, 111_236_476, "ISCA-37405", "2G-pos-1")] // adds LIMS4
    fn handle_cases_2d_2g(
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

        let res = evaluator.handle_cases_2d_2g(&sv_interval, &benign_regions)?;
        insta::assert_yaml_snapshot!(res);

        Ok(())
    }

    /// Test internal working of `handle_cases_2h_l`.
    ///
    /// We cannot match phenotypes, so we never return 2K and just test 2J.
    #[tracing_test::traced_test]
    #[rstest::rstest]
    // Note: same cases as in `evaluate` above -- keep in sync!
    //
    // Case 2H: HI gene fully contained within observed copy number gain.
    #[case("2", 189839099, 189877472, "COL3A1", "2H-pos-1")]
    // Case 2I: Both breakpoints in the same HI gene.
    #[case("2", 189839100, 189877471, "COL3A1", "2I-pos-1")]
    // Case 2J: One breakpoint in an HI gene.
    #[case("2", 189877072, 189877500, "COL3A1", "2J-pos-1")] // left bp
    #[case("2", 189839000, 189839599, "COL3A1", "2J-pos-2")]
    // right bp
    // No Case 2K: we cannot handle phenotypes well.
    // Case 2L: One breakpoint in any gene.
    #[case("6", 111_620_000, 111_621_000, "REV3L", "2L-pos-1")] // REV3L left bp
    #[case("6", 111_800_000, 111_805_000, "REV3L", "2L-pos-1")] // REV3L right bp
    fn handle_cases_2h_2l(
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

        let (clingen_genes, _) = global_evaluator_37.clingen_overlaps(&sv_interval);
        let res = evaluator.handle_cases_2h_2l(&sv_interval, &clingen_genes)?;
        insta::assert_yaml_snapshot!(res);

        Ok(())
    }
}
