//! Data structures for representing the actual results of DEL/copy number loss.

use strum::IntoEnumIterator as _;

use crate::strucvars::{
    data::{clingen_dosage, hgnc::GeneIdInfo},
    eval::{
        common::{FunctionalElement, GeneOverlap, HasScoreRange, SuggestedScore},
        result::{ClinicalSignificance, ClinvarSvOverlap, GnomadSvOverlap, Pvs1Result, ScoreRange},
    },
};

/// Evaluation summary.
#[derive(Debug, Clone, Default, PartialEq, serde::Deserialize, serde::Serialize)]
pub struct Summary {
    /// Overall score.
    pub score: f32,
    /// Overall clinical significance.
    pub clinical_significance: ClinicalSignificance,
}

impl Summary {
    /// Construct with sections.
    pub fn with_sections(sections: &Vec<Section>) -> Self {
        let score = Self::get_score(sections);
        Self {
            score,
            clinical_significance: Self::get_clinical_signifiance(score),
        }
    }

    /// Return overall score.
    fn get_score(sections: &Vec<Section>) -> f32 {
        sections.iter().map(|s| s.suggested_score()).sum()
    }

    /// Return overall clinical significance.
    fn get_clinical_signifiance(score: f32) -> ClinicalSignificance {
        for clinsig in ClinicalSignificance::iter() {
            let itv = clinsig.range();
            if itv.contains(&ScoreRange::Singleton { at: score }) {
                return clinsig;
            }
        }

        unreachable!()
    }
}

/// Overall evaluation results for all sections.
#[derive(Debug, Clone, Default, PartialEq, serde::Deserialize, serde::Serialize)]
pub struct Evaluation {
    /// Summary of the evaluation.
    pub summary: Summary,
    /// All section results.
    pub sections: Vec<Section>,
}

impl Evaluation {
    /// Construct with sections.
    pub fn from_sections(sections: Vec<Section>) -> Self {
        let summary = Summary::with_sections(&sections);
        Self { sections, summary }
    }
}

impl Evaluation {
    /// Return overall score.
    pub fn score(&self) -> f32 {
        self.sections.iter().map(|s| s.suggested_score()).sum()
    }

    /// Return the overall clinical significance.
    pub fn clinical_signifiance(&self) -> ClinicalSignificance {
        let score = self.score();
        for clinsig in ClinicalSignificance::iter() {
            let itv = clinsig.range();
            if itv.contains(&ScoreRange::Singleton { at: score }) {
                return clinsig;
            }
        }

        unreachable!()
    }
}

/// Evaluation results for each section of the ACMG rule.
#[derive(Debug, Clone, PartialEq, serde::Deserialize, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Section {
    /// Results of Section L1.
    L1(L1),
    /// Results of Section L2.
    L2(L2),
    /// Results of Section L3.
    L3(L3),
    /// Results of Section L4.
    L4(L4),
}

impl SuggestedScore for Section {
    fn suggested_score(&self) -> f32 {
        match self {
            Section::L1(L1::L1A(l1a)) => l1a.suggested_score(),
            Section::L1(L1::L1B(l1b)) => l1b.suggested_score(),
            Section::L2(L2::L2A(l2a)) => l2a.suggested_score(),
            Section::L2(L2::L2B(l2b)) => l2b.suggested_score(),
            Section::L2(L2::L2C1(l2c1)) => l2c1.suggested_score(),
            Section::L2(L2::L2C2(l2c2)) => l2c2.suggested_score(),
            Section::L2(L2::L2D1(l2d1)) => l2d1.suggested_score(),
            Section::L2(L2::L2D2(l2d2)) => l2d2.suggested_score(),
            Section::L2(L2::L2D3(l2d3)) => l2d3.suggested_score(),
            Section::L2(L2::L2D4(l2d4)) => l2d4.suggested_score(),
            Section::L2(L2::L2E(l2e)) => l2e.suggested_score(),
            Section::L2(L2::L2F(l2f)) => l2f.suggested_score(),
            Section::L2(L2::L2G(l2g)) => l2g.suggested_score(),
            Section::L2(L2::L2H(l2h)) => l2h.suggested_score(),
            Section::L3(L3::L3A(l3a)) => l3a.suggested_score(),
            Section::L3(L3::L3B(l3b)) => l3b.suggested_score(),
            Section::L3(L3::L3C(l3c)) => l3c.suggested_score(),
            Section::L4(L4::L4Patho(l4p)) => l4p.suggested_score(),
            Section::L4(L4::L4N(l4b)) => l4b.suggested_score(),
            Section::L4(L4::L4O(l4o)) => l4o.suggested_score(),
        }
    }
}

impl HasScoreRange for Section {
    fn min_score(&self) -> f32 {
        match self {
            Section::L1(L1::L1A(_)) => 0.0,
            Section::L1(L1::L1B(_)) => -0.6,
            Section::L2(L2::L2A(_)) => 1.0,
            Section::L2(L2::L2B(_)) => 0.0,
            Section::L2(L2::L2C1(_)) => 0.9,
            Section::L2(L2::L2C2(_)) => 0.0,
            Section::L2(L2::L2D1(_)) => 0.0,
            Section::L2(L2::L2D2(_)) => 0.45,
            Section::L2(L2::L2D3(_)) => 0.0,
            Section::L2(L2::L2D4(_)) => 0.45,
            Section::L2(L2::L2E(_)) => 0.0,
            Section::L2(L2::L2F(_)) => -1.0,
            Section::L2(L2::L2G(_)) => 0.0,
            Section::L2(L2::L2H(_)) => 0.15,
            Section::L3(L3::L3A(_)) => 0.0,
            Section::L3(L3::L3B(_)) => 0.45,
            Section::L3(L3::L3C(_)) => 0.9,
            Section::L4(L4::L4Patho(_)) => -0.9,
            Section::L4(L4::L4N(_)) => -0.9,
            Section::L4(L4::L4O(_)) => -1.0,
        }
    }

    fn max_score(&self) -> f32 {
        match self {
            Section::L1(L1::L1A(_)) => 0.0,
            Section::L1(L1::L1B(_)) => -0.6,
            Section::L2(L2::L2A(_)) => 1.0,
            Section::L2(L2::L2B(_)) => 0.0,
            Section::L2(L2::L2C1(_)) => 1.0,
            Section::L2(L2::L2C2(_)) => 0.45,
            Section::L2(L2::L2D1(_)) => 0.0,
            Section::L2(L2::L2D2(_)) => 0.90,
            Section::L2(L2::L2D3(_)) => 0.45,
            Section::L2(L2::L2D4(_)) => 1.0,
            Section::L2(L2::L2E(_)) => 0.9,
            Section::L2(L2::L2F(_)) => -1.0,
            Section::L2(L2::L2G(_)) => 0.0,
            Section::L2(L2::L2H(_)) => 0.15,
            Section::L3(L3::L3A(_)) => 0.0,
            Section::L3(L3::L3B(_)) => 0.45,
            Section::L3(L3::L3C(_)) => 0.9,
            Section::L4(L4::L4Patho(_)) => 0.9,
            Section::L4(L4::L4N(_)) => 0.0,
            Section::L4(L4::L4O(_)) => 0.0,
        }
    }
}

/// Enumeration of the categories for the copy number loss evaluation, Section 1.
#[derive(Debug, Clone, PartialEq, serde::Deserialize, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum L1 {
    /// Contains protein-coding or other known functionally important elements.
    L1A(L1A),
    /// Does NOT contain protein-coding or any functionally important elements.
    L1B(L1B),
}

/// Result of the L1A subsection (important feature).
#[derive(Debug, Clone, PartialEq, Default, serde::Deserialize, serde::Serialize)]
pub struct L1A {
    /// Overlapping transcripts/genes.
    pub genes: Vec<GeneOverlap>,
    /// Overlapping functional elements.
    pub functional_elements: Vec<FunctionalElement>,
}

impl SuggestedScore for L1A {
    fn suggested_score(&self) -> f32 {
        0.0
    }
}

/// Result of the L1B subsection (no important feature).
#[derive(Debug, Clone, PartialEq, Default, serde::Deserialize, serde::Serialize)]
pub struct L1B {
    // no members
}

impl SuggestedScore for L1B {
    fn suggested_score(&self) -> f32 {
        -0.6
    }
}

/// Enumeration of the categories for the structural variant evaluation, Section 2.
#[derive(Debug, Clone, PartialEq, serde::Deserialize, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum L2 {
    /// Complete overlap of an established HI gene/genomic region.
    L2A(L2A),
    /// Partial overlap of an established HI genomic region.
    L2B(L2B),
    /// Partial overlap of 5' end of established HI gene ... coding sequence is involved.
    L2C1(L2C1),
    /// Partial overlap of 5' end of established HI gene ... only the 5' UTR is involved.
    L2C2(L2C2),
    /// Partial overlap of 3' end of established HI gene ... only the 3' UTR is involved.
    L2D1(L2D1),
    /// Partial overlap of 3' end of established HI gene ... only the last exon is involved,
    /// etablished pathogenic variants exist in this exon.
    L2D2(L2D2),
    /// Partial overlap of 3' end of established HI gene ... only the last exon is involved,
    /// no established pathogenic variants exist in this exon.
    L2D3(L2D3),
    /// Partial overlap of 3' end of established HI gene ... it includes other exons in addition
    /// to the last exon; nonsense-mediated decay is expected to occur.
    L2D4(L2D4),
    /// Both breakpoints are within the same gene (intragenic CNV, gene-level sequence variant
    /// ... PVS1 rules apply).
    L2E(L2E),
    /// Completely contained within an established benign CNV region.
    L2F(L2F),
    /// Overlaps with an established benign CNV, but includes additional genomic material.
    L2G(L2G),
    /// Two or more HI predictors suggest that AT LEAST ONE gene in the interval is HI.
    L2H(L2H),
}

/// Result of the L2A subsection (complete overlap of an established HI gene/genomic region).
#[derive(Debug, Clone, PartialEq, Default, serde::Deserialize, serde::Serialize)]
pub struct L2A {
    /// Suggested score for the subsection.
    pub suggested_score: f32,
    /// Overlapping HI genes.
    pub hi_genes: Vec<clingen_dosage::Gene>,
    /// Overlapping HI genomic regions.
    pub hi_regions: Vec<clingen_dosage::Region>,
}

impl SuggestedScore for L2A {
    fn suggested_score(&self) -> f32 {
        self.suggested_score
    }
}

/// Result of the L2B subsection (partial overlap of an established HI genomic region).
#[derive(Debug, Clone, PartialEq, Default, serde::Deserialize, serde::Serialize)]
pub struct L2B {
    // no members, just signal to the user
}

impl SuggestedScore for L2B {
    fn suggested_score(&self) -> f32 {
        0.0
    }
}

/// Result of the L2C1 subsection (partial overlap of 5' end of established HI gene
/// ... coding sequence is involved).
#[derive(Debug, Clone, PartialEq, serde::Deserialize, serde::Serialize)]
pub struct L2C1 {
    /// Suggested score for the subsection.
    pub suggested_score: f32,
    /// Overlapping gene information.
    pub gene: GeneIdInfo,
}

impl SuggestedScore for L2C1 {
    fn suggested_score(&self) -> f32 {
        self.suggested_score
    }
}

/// Result of the L2C2 subsection (partial overlap of 5' end of established HI gene
/// ... only the 5' UTR is involved).
#[derive(Debug, Clone, PartialEq, serde::Deserialize, serde::Serialize)]
pub struct L2C2 {
    /// Suggested score for the subsection.
    pub suggested_score: f32,
    /// Overlapping gene information.
    pub gene: GeneIdInfo,
}

impl SuggestedScore for L2C2 {
    fn suggested_score(&self) -> f32 {
        self.suggested_score
    }
}

/// Result of the L2D1 subsection (partial overlap of 3' end of established HI gene
/// ... only the 3' UTR is involved).
#[derive(Debug, Clone, PartialEq, serde::Deserialize, serde::Serialize)]
pub struct L2D1 {
    /// Overlapping gene information.
    pub gene: GeneIdInfo,
}

impl SuggestedScore for L2D1 {
    fn suggested_score(&self) -> f32 {
        0.0
    }
}

/// Result of the L2D2 subsection (partial overlap of 3' end of established HI gene
/// ... only the last exon is involved; other established pathogenic variants have
/// been reported in this exon).
#[derive(Debug, Clone, PartialEq, serde::Deserialize, serde::Serialize)]
pub struct L2D2 {
    /// Suggested score for the subsection.
    pub suggested_score: f32,
    /// Overlappign gene information.
    pub gene: GeneIdInfo,
    /// VCV identifiers of ClinVar variants supporting the score.
    pub clinvar_ids: Vec<String>,
}

impl SuggestedScore for L2D2 {
    fn suggested_score(&self) -> f32 {
        self.suggested_score
    }
}

/// Result of the L2D3 subsection (partial overlap of 3' end of established HI gene
/// ... only the last exon is involved; no other established pathogenic variants have
///  been reported in this exon).
#[derive(Debug, Clone, PartialEq, Default, serde::Deserialize, serde::Serialize)]
pub struct L2D3 {
    /// Suggested score for the subsection.
    pub suggested_score: f32,
    /// Overlapping gene information.
    pub gene: GeneIdInfo,
}

impl SuggestedScore for L2D3 {
    fn suggested_score(&self) -> f32 {
        self.suggested_score
    }
}

/// Result of the L2D4 subsection (partial overlap of 3' end of established HI gene
/// ... it includes other exons in addition to the last exon; nonsense-mediate
/// decay is expected to occur).
#[derive(Debug, Clone, PartialEq, Default, serde::Deserialize, serde::Serialize)]
pub struct L2D4 {
    /// Suggested score for the subsection.
    pub suggested_score: f32,
    /// Overlapping gene information.
    pub gene: GeneIdInfo,
    /// Numbers of the overlapping exons.
    pub exon_nos: Vec<u32>,
}

impl SuggestedScore for L2D4 {
    fn suggested_score(&self) -> f32 {
        self.suggested_score
    }
}

/// Result of the L2E subsection (both breakpoints are within the same gene (intragenic
/// CNV, gene-level sequence variant). PVS1 rules apply from ClinGen SVI WG.
#[derive(Debug, Clone, PartialEq, Default, serde::Deserialize, serde::Serialize)]
pub struct L2E {
    /// Suggested score for the subsection.
    pub suggested_score: f32,
    /// PVS1 assessment result.
    pub pvs1_result: Pvs1Result,
}

impl SuggestedScore for L2E {
    fn suggested_score(&self) -> f32 {
        self.suggested_score
    }
}

/// Result of the L2F subsection (completely contained with an established benign CNV region).
#[derive(Debug, Clone, PartialEq, Default, serde::Deserialize, serde::Serialize)]
pub struct L2F {
    /// Suggested score for the subsection.
    pub suggested_score: f32,
    /// Overlapping CNV benign genomic regions.
    pub benign_regions: Vec<clingen_dosage::Region>,
}

impl SuggestedScore for L2F {
    fn suggested_score(&self) -> f32 {
        self.suggested_score
    }
}

/// Result of the L2G subsection (overlaps with an established benign CNV, but includes additional
///  genomic material.
#[derive(Debug, Clone, PartialEq, Default, serde::Deserialize, serde::Serialize)]
pub struct L2G {
    /// Suggested score for the subsection.
    pub suggested_score: f32,
    /// Overlapping CNV benign genomic regions.
    pub benign_regions: Vec<clingen_dosage::Region>,
}

impl SuggestedScore for L2G {
    fn suggested_score(&self) -> f32 {
        self.suggested_score
    }
}

/// HI predictor results for use in `L2H`.
#[derive(Debug, Clone, PartialEq, serde::Deserialize, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum GeneHiPredictorResult {
    /// gnomAD pLI results
    GnomadPli {
        /// gnomAD pLI score
        pli_score: f64,
        /// upper bound of LoF observed/expected confidence interval
        oe_lof_upper: f64,
    },
    /// DECIPHER HI index
    DecipherHiIndex {
        /// DECIPHER HI index
        hi_index: f64,
    },
}

/// HI prediction results for a single gene.
#[derive(Debug, Clone, PartialEq, serde::Deserialize, serde::Serialize)]
pub struct GeneHiPrediction {
    /// Identifier of the gene.
    pub gene: GeneIdInfo,
    /// Results of the HI predictors.
    pub results: Vec<GeneHiPredictorResult>,
}

/// Result of the L2H subsection (two or more HI predictors suggest that AT LEAST ONE gene in the
/// interval is HI).
#[derive(Debug, Clone, PartialEq, Default, serde::Deserialize, serde::Serialize)]
pub struct L2H {
    /// Suggested score for the subsection.
    pub suggested_score: f32,
    /// HGNC identifiers of the overlapping predicted HI genes.
    pub gene_hi_predictions: Vec<GeneHiPrediction>,
}

impl SuggestedScore for L2H {
    fn suggested_score(&self) -> f32 {
        self.suggested_score
    }
}

/// Enumeration of the categories for the structural variant evaluation, Section 2.
#[derive(Debug, Clone, PartialEq, serde::Deserialize, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum L3 {
    /// <=24 genes.
    L3A(L3Count),
    /// 25-34 genes.
    L3B(L3Count),
    /// >=35 genes.
    L3C(L3Count),
}

/// Result of the 3A subsection (Number of protein-coding RefSeq genes wholly or partially included
/// in the copy-number loss).
#[derive(Debug, Clone, PartialEq, Default, serde::Deserialize, serde::Serialize)]
pub struct L3Count {
    /// Suggested score for the subsection.
    pub suggested_score: f32,
    /// Number of protein-coding RefSeq genes wholly or partially included in the copy-number loss.
    pub num_genes: usize,
    /// The overlapping genes.
    pub genes: Vec<GeneOverlap>,
}

impl SuggestedScore for L3Count {
    fn suggested_score(&self) -> f32 {
        self.suggested_score
    }
}

/// Enumeration of the categories for the structural variant evaluation, Section 4.
///
/// Only 4O can be automatically determined.
#[derive(Debug, Clone, PartialEq, serde::Deserialize, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum L4 {
    /// Overlap with pathogenic variants; must be evaluated by a human.
    ///
    /// This could be one of 4A, 4B, 4C, 4D, 4E, 4L, 4M.
    L4Patho(L4Patho),
    /// Overlap with benign variants; must be evaluated by a human.
    ///
    /// This is roughly equiavalent to 4N.
    L4N(L4N),
    /// Case–control and population evidence; Overlap with common population variation.
    L4O(L4O),
}

/// Pathogenic variants from ClinVar, must be resolved by a human.
///
/// This could be 4D, 4E, 4M.
#[derive(Debug, Clone, PartialEq, Default, serde::Deserialize, serde::Serialize)]
pub struct L4Patho {
    /// Overlapping ClinVar variants, descendingly by overlap.
    pub overlaps: Vec<ClinvarSvOverlap>,
}

impl SuggestedScore for L4Patho {
    fn suggested_score(&self) -> f32 {
        0.0
    }
}

/// Benign variants from ClinVar.
///
/// This is roughly 4N.
#[derive(Debug, Clone, PartialEq, Default, serde::Deserialize, serde::Serialize)]
pub struct L4N {
    /// Overlapping ClinVar variants, descendingly by overlap.
    pub overlaps: Vec<ClinvarSvOverlap>,
}

impl SuggestedScore for L4N {
    fn suggested_score(&self) -> f32 {
        -0.9
    }
}

/// Result of the 4O subsection (overlap with common population variation).
#[derive(Debug, Clone, PartialEq, Default, serde::Deserialize, serde::Serialize)]
pub struct L4O {
    /// Suggested score for the subsection.
    pub suggested_score: f32,
    /// Accession identifiers of overlapping common variants.
    pub overlaps: Vec<GnomadSvOverlap>,
}

impl SuggestedScore for L4O {
    fn suggested_score(&self) -> f32 {
        self.suggested_score
    }
}
