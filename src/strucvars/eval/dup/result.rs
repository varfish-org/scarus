//! Data structures for representing the actual results of DEL/copy number gain.

use crate::strucvars::{
    data::{clingen_dosage, hgnc::GeneIdInfo},
    eval::{
        common::{FunctionalElement, GeneOverlap, ScoreRange, SuggestedScore},
        result::{ClinvarSvOverlap, GnomadSvOverlap, Pvs1Result},
    },
};

/// Evaluation results for each section of the ACMG rule.
#[derive(Debug, Clone, PartialEq, serde::Deserialize, serde::Serialize)]
pub enum Section {
    /// Results of Section G1.
    G1(G1),
    /// Results of Section G2.
    G2(G2),
    /// Results of Section G3.
    G3(G3),
    /// Results of Section G4.
    G4(G4),
}

impl SuggestedScore for Section {
    fn suggested_score(&self) -> f32 {
        match self {
            Section::G1(G1::G1A(g1a)) => g1a.suggested_score(),
            Section::G1(G1::G1B(g1b)) => g1b.suggested_score(),
            Section::G2(G2::G2A(g2a)) => g2a.suggested_score(),
            Section::G2(G2::G2B(g2b)) => g2b.suggested_score(),
            Section::G2(G2::G2C(g2c)) => g2c.suggested_score(),
            Section::G2(G2::G2D(g2d)) => g2d.suggested_score(),
            Section::G2(G2::G2E(g2e)) => g2e.suggested_score(),
            Section::G2(G2::G2F(g2f)) => g2f.suggested_score(),
            Section::G2(G2::G2G(g2g)) => g2g.suggested_score(),
            Section::G2(G2::G2H(g2h)) => g2h.suggested_score(),
            Section::G2(G2::G2I(g2i)) => g2i.suggested_score(),
            Section::G2(G2::G2J(g2j)) => g2j.suggested_score(),
            Section::G2(G2::G2L(g2l)) => g2l.suggested_score(),
            Section::G3(G3::G3A(g3a)) => g3a.suggested_score(),
            Section::G3(G3::G3B(g3b)) => g3b.suggested_score(),
            Section::G3(G3::G3C(g3c)) => g3c.suggested_score(),
            Section::G4(G4::G4Patho(g4p)) => g4p.suggested_score(),
            Section::G4(G4::G4N(g4n)) => g4n.suggested_score(),
            Section::G4(G4::G4O(g4o)) => g4o.suggested_score(),
        }
    }
}

impl ScoreRange for Section {
    fn min_score(&self) -> f32 {
        match self {
            Section::G1(G1::G1A(_)) => 0.0,
            Section::G1(G1::G1B(_)) => -0.6,
            Section::G2(G2::G2A(_)) => 1.0,
            Section::G2(G2::G2B(_)) => 0.0,
            Section::G2(G2::G2C(_)) => -1.0,
            Section::G2(G2::G2D(_)) => -1.0,
            Section::G2(G2::G2E(_)) => 0.0,
            Section::G2(G2::G2F(_)) => -1.0,
            Section::G2(G2::G2G(_)) => 0.0,
            Section::G2(G2::G2H(_)) => 0.0,
            Section::G2(G2::G2I(_)) => 0.0,
            Section::G2(G2::G2J(_)) => 0.0,
            Section::G2(G2::G2L(_)) => 0.0,
            Section::G3(G3::G3A(_)) => 0.0,
            Section::G3(G3::G3B(_)) => 0.45,
            Section::G3(G3::G3C(_)) => 0.9,
            Section::G4(G4::G4Patho(_)) => -0.9,
            Section::G4(G4::G4N(_)) => -0.9,
            Section::G4(G4::G4O(_)) => -1.0,
        }
    }

    fn max_score(&self) -> f32 {
        match self {
            Section::G1(G1::G1A(_)) => 0.0,
            Section::G1(G1::G1B(_)) => -0.6,
            Section::G2(G2::G2A(_)) => 1.0,
            Section::G2(G2::G2B(_)) => 0.0,
            Section::G2(G2::G2C(_)) => -1.0,
            Section::G2(G2::G2D(_)) => -1.0,
            Section::G2(G2::G2E(_)) => 0.0,
            Section::G2(G2::G2F(_)) => 0.0,
            Section::G2(G2::G2G(_)) => 0.0,
            Section::G2(G2::G2H(_)) => 0.0,
            Section::G2(G2::G2I(_)) => 0.9,
            Section::G2(G2::G2J(_)) => 0.0,
            Section::G2(G2::G2L(_)) => 0.0,
            Section::G3(G3::G3A(_)) => 0.0,
            Section::G3(G3::G3B(_)) => 0.45,
            Section::G3(G3::G3C(_)) => 0.9,
            Section::G4(G4::G4Patho(_)) => 0.9,
            Section::G4(G4::G4N(_)) => 0.9,
            Section::G4(G4::G4O(_)) => 0.0,
        }
    }
}

/// Enumeration of the categories for the copy number gain evaluation, Section 1.
#[derive(Debug, Clone, PartialEq, serde::Deserialize, serde::Serialize)]
pub enum G1 {
    /// Contains protein-coding or other known functionally important elements.
    G1A(G1A),
    /// Does NOT contain protein-coding or any functionally important elements.
    G1B(G1B),
}

/// Result of the G1A subsection (important feature).
#[derive(Debug, Clone, PartialEq, Default, serde::Deserialize, serde::Serialize)]
pub struct G1A {
    /// Overlapping transcripts/genes.
    pub genes: Vec<GeneOverlap>,
    /// Overlapping functional elements.
    pub functional_elements: Vec<FunctionalElement>,
}

impl SuggestedScore for G1A {
    fn suggested_score(&self) -> f32 {
        0.0
    }
}

/// Result of the G1B subsection (no important feature).
#[derive(Debug, Clone, PartialEq, Default, serde::Deserialize, serde::Serialize)]
pub struct G1B {
    // no members
}

impl SuggestedScore for G1B {
    fn suggested_score(&self) -> f32 {
        -0.6
    }
}

/// Enumeration of the categories for the structural variant evaluation, Section 2.
#[derive(Debug, Clone, PartialEq, serde::Deserialize, serde::Serialize)]
pub enum G2 {
    /// Complete overlap; Tthe TS gene or minimal criticala reegion is fully contained
    /// within the observed copy number gain.
    G2A(G2A),
    /// Partial overlap of an established TS region.
    G2B(G2B),
    /// Identical in gene content to the established benign copy number gain.
    G2C(G2C),
    /// Smaller than established benign copy number gain, breakpoint(s) does not
    /// interrupt protein-coding genes.
    G2D(G2D),
    /// Smaller than established benign copy number gain, breakpoint(s) potentially
    /// interrupt protein-coding genes.
    G2E(G2E),
    /// Larger tahn known benign copy-number gain, does not include additional
    /// protein-coding genes.
    G2F(G2F),
    /// Overlaps a benign copy number gain but includes additional material.
    G2G(G2G),
    /// TS gene fully contained within obsrved copy-number gain.
    G2H(G2H),
    /// Both breakpoints are within the asme gene (gene-level sequence variant, possibly
    /// resulting in LoF.
    G2I(G2I),
    /// One breakpoint is within an established TS gene, patient's phenotype is either
    /// inconsistent with what is expecte for LOF or that gene OR unknown
    G2J(G2J),
    /// One or both breakpoints are within gene(s) of no established clinical significance.
    G2L(G2L),
}

/// Result of the G2A subsection (complete overlap).
#[derive(Debug, Clone, PartialEq, Default, serde::Deserialize, serde::Serialize)]
pub struct G2A {
    /// Suggested score for the subsection.
    pub suggested_score: f32,
    /// Overlapping TS genes.
    pub ts_genes: Vec<clingen_dosage::Gene>,
    /// Overlapping TS genomic regions.
    pub ts_regions: Vec<clingen_dosage::Region>,
}

impl SuggestedScore for G2A {
    fn suggested_score(&self) -> f32 {
        self.suggested_score
    }
}

/// Result of the G2B subsection (partial overlap of established TS region).
#[derive(Debug, Clone, PartialEq, Default, serde::Deserialize, serde::Serialize)]
pub struct G2B {
    /// Overlapping TS genomic regions.
    pub ts_regions: Vec<clingen_dosage::Region>,
}

impl SuggestedScore for G2B {
    fn suggested_score(&self) -> f32 {
        0.0
    }
}

/// Result of the G2C subsection (identical in gene content to the established
/// benign TS region).
#[derive(Debug, Clone, PartialEq, serde::Deserialize, serde::Serialize)]
pub struct G2C {
    /// Suggested score for the subsection.
    pub suggested_score: f32,
    /// Overlapping TS genomic regions.
    pub benign_regions: Vec<clingen_dosage::Region>,
}

impl SuggestedScore for G2C {
    fn suggested_score(&self) -> f32 {
        self.suggested_score
    }
}

/// Result of the G2D subsection (smaller than established benign copy-number gain,
/// breakpoint(s) does not interrupt protein-coding gene).
#[derive(Debug, Clone, PartialEq, serde::Deserialize, serde::Serialize)]
pub struct G2D {
    /// Suggested score for the subsection.
    pub suggested_score: f32,
    /// Overlapping benign genomic region.
    pub benign_region: clingen_dosage::Region,
}

impl SuggestedScore for G2D {
    fn suggested_score(&self) -> f32 {
        self.suggested_score
    }
}

/// Result of the G2E subsection (smaller than established benign copy-number gain,
/// breakpoint(s) potentially interrupt protein-coding gene).
#[derive(Debug, Clone, PartialEq, serde::Deserialize, serde::Serialize)]
pub struct G2E {
    /// Suggested score for the subsection.
    pub suggested_score: f32,
    /// Overlapping benign genomic regions.
    pub benign_region: clingen_dosage::Region,
    /// Potentially affected genes.
    pub genes: Vec<GeneIdInfo>,
    /// Overlapping functional elements.
    pub functional_elements: Vec<FunctionalElement>,
}

impl SuggestedScore for G2E {
    fn suggested_score(&self) -> f32 {
        self.suggested_score
    }
}

/// Result of G2F subsection (larger than known benign copy-number gain,
/// does not include additional protein-coding genes).
#[derive(Debug, Clone, PartialEq, serde::Deserialize, serde::Serialize)]
pub struct G2F {
    /// Suggested score for the subsection.
    pub suggested_score: f32,
    /// Overlapping benign genomic regions.
    pub benign_region: clingen_dosage::Region,
}

impl SuggestedScore for G2F {
    fn suggested_score(&self) -> f32 {
        self.suggested_score
    }
}

/// Result of the G2G subsection (overlaps a benign copy-number gain but includes
/// additional genomic material).
#[derive(Debug, Clone, PartialEq, serde::Deserialize, serde::Serialize)]
pub struct G2G {
    /// Suggested score for the subsection.
    pub suggested_score: f32,
    /// Overlapping benign genomic regions.
    pub benign_region: clingen_dosage::Region,
    /// Additional genes.
    pub genes: Vec<GeneIdInfo>,
    /// Overlapping functional elements.
    pub functional_elements: Vec<FunctionalElement>,
}

impl SuggestedScore for G2G {
    fn suggested_score(&self) -> f32 {
        self.suggested_score
    }
}

/// Result of the G2H subsection (HI gene fully contained within observed copy-number
/// gain).
#[derive(Debug, Clone, PartialEq, Default, serde::Deserialize, serde::Serialize)]
pub struct G2H {
    /// Suggested score for the subsection.
    pub suggested_score: f32,
    /// Overlapping HI gene information.
    pub hi_genes: Vec<GeneIdInfo>,
}

impl SuggestedScore for G2H {
    fn suggested_score(&self) -> f32 {
        self.suggested_score
    }
}

/// Result of the G2I subsection (both breakpoint are within the same gene, potential
/// LOF, fall back to PVS1).
#[derive(Debug, Clone, PartialEq, Default, serde::Deserialize, serde::Serialize)]
pub struct G2I {
    /// Suggested score for the subsection.
    pub suggested_score: f32,
    /// PVS assessment result.
    pub pvs1_result: Pvs1Result,
    /// HI genes for which the assessment was made.
    pub hi_genes: Vec<GeneIdInfo>,
}

impl SuggestedScore for G2I {
    fn suggested_score(&self) -> f32 {
        self.suggested_score
    }
}

/// Result of the G2J subsection (One breakpoint is within an esbalished HI gene,
/// patient's phenotype is either inconsistent with what is expected for LOF of
/// that gene OR unknown).
#[derive(Debug, Clone, PartialEq, Default, serde::Deserialize, serde::Serialize)]
pub struct G2J {
    /// Suggested score for the subsection.
    pub suggested_score: f32,
    /// HI genes overlapping with the CNV.
    pub hi_genes: Vec<GeneIdInfo>,
}

impl SuggestedScore for G2J {
    fn suggested_score(&self) -> f32 {
        self.suggested_score
    }
}

/// Result of the G2L subsection (one or both breakpoints are within gene(s) of
/// no established linical significance).
#[derive(Debug, Clone, PartialEq, Default, serde::Deserialize, serde::Serialize)]
pub struct G2L {
    /// Suggested score for the subsection.
    pub suggested_score: f32,
    /// Genes overlapping with the DUP.
    pub genes: Vec<GeneIdInfo>,
}

impl SuggestedScore for G2L {
    fn suggested_score(&self) -> f32 {
        self.suggested_score
    }
}

/// Enumeration of the categories for the structural variant evaluation, Section 2.
#[derive(Debug, Clone, PartialEq, serde::Deserialize, serde::Serialize)]
pub enum G3 {
    /// <=34 genes.
    G3A(G3Count),
    /// 35-49 genes.
    G3B(G3Count),
    /// >=50 genes.
    G3C(G3Count),
}

/// Result of the 3A subsection (Number of protein-coding RefSeq genes wholly or partially included
/// in the copy-number loss).
#[derive(Debug, Clone, PartialEq, Default, serde::Deserialize, serde::Serialize)]
pub struct G3Count {
    /// Suggested score for the subsection.
    pub suggested_score: f32,
    /// Number of protein-coding RefSeq genes wholly or partially included in the copy-number loss.
    pub num_genes: usize,
    /// The overlapping genes.
    pub genes: Vec<GeneOverlap>,
}

impl SuggestedScore for G3Count {
    fn suggested_score(&self) -> f32 {
        self.suggested_score
    }
}

/// Enumeration of the categories for the structural variant evaluation, Section 4.
///
/// Only 4O can be automatically determined.
#[derive(Debug, Clone, PartialEq, serde::Deserialize, serde::Serialize)]
pub enum G4 {
    /// Overlap with pathogenic variants; must be evaluated by a human.
    ///
    /// This could be one of 4A, 4B, 4C, 4D, 4E, 4L, 4M.
    G4Patho(G4Patho),
    /// Overlap with benign variants; must be evaluated by a human.
    ///
    /// This is roughly equiavalent to 4N.
    G4N(G4N),
    /// Case–control and population evidence; Overlap with common population variation.
    G4O(G4O),
}

/// Result of the 4O subsection (overlap with common population variation).
#[derive(Debug, Clone, PartialEq, Default, serde::Deserialize, serde::Serialize)]
pub struct G4Patho {
    /// Overlapping ClinVar variants, descendingly by overlap.
    pub overlaps: Vec<ClinvarSvOverlap>,
}

impl SuggestedScore for G4Patho {
    fn suggested_score(&self) -> f32 {
        0.0
    }
}

/// Benign variants from ClinVar.
///
/// This is roughly 4N.
#[derive(Debug, Clone, PartialEq, Default, serde::Deserialize, serde::Serialize)]
pub struct G4N {
    /// Overlapping ClinVar variants, descendingly by overlap.
    pub overlaps: Vec<ClinvarSvOverlap>,
}

impl SuggestedScore for G4N {
    fn suggested_score(&self) -> f32 {
        -0.9
    }
}

/// Result of the 4O subsection (overlap with common population variation).
#[derive(Debug, Clone, PartialEq, Default, serde::Deserialize, serde::Serialize)]
pub struct G4O {
    /// Suggested score for the subsection.
    pub suggested_score: f32,
    /// Accession identifiers of overlapping common variants.
    pub overlaps: Vec<GnomadSvOverlap>,
}

impl SuggestedScore for G4O {
    fn suggested_score(&self) -> f32 {
        self.suggested_score
    }
}
