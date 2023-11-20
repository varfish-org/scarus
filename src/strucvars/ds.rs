//! Shared data structures for `strucvars`.

use super::data::intervals::Interval;

/// Enumeration for SV type.
#[derive(
    Debug,
    Clone,
    Copy,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Hash,
    Default,
    serde::Deserialize,
    serde::Serialize,
)]
#[serde(rename_all = "snake_case")]
pub enum SvType {
    /// Deletion, copy number loss.
    #[default]
    Del,
    /// Duplication, copy number gain.
    Dup,
}

/// Representation of inner / outer coordinates.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, serde::Deserialize, serde::Serialize)]
pub struct AmbiguousRange {
    /// Outer start position.
    pub outer_start: u32,
    /// Inner start position.
    pub inner_start: u32,
    /// Inner end position.
    pub inner_end: u32,
    /// Outer end position.
    pub outer_end: u32,
}

/// Representation of a structural variant.
#[derive(Debug, Clone, PartialEq, Default, serde::Deserialize, serde::Serialize)]
pub struct StructuralVariant {
    /// Chromosome.
    pub chrom: String,
    /// 1-based start position.
    pub start: u32,
    /// 1-based stop position.
    pub stop: u32,
    /// SV type.
    pub svtype: SvType,

    /// Optional ambiguous start/end positions.
    pub ambiguous_range: Option<AmbiguousRange>,
}

impl From<StructuralVariant> for Interval {
    fn from(val: StructuralVariant) -> Self {
        Interval::new(
            val.chrom,
            val.start.saturating_sub(1).into()..val.stop.into(),
        )
    }
}
