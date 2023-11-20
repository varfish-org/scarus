//! Shared data structures for `strucvars`.

use std::{fmt::Display, str::FromStr};

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

impl Display for SvType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SvType::Del => write!(f, "DEL"),
            SvType::Dup => write!(f, "DUP"),
        }
    }
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

impl FromStr for StructuralVariant {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let arr = s.split(':').collect::<Vec<_>>();
        if arr.len() != 4 {
            anyhow::bail!("invalid number of fields: {:?}", &arr);
        }
        let svtype = match arr[0] {
            "DEL" => SvType::Del,
            "DUP" => SvType::Dup,
            _ => anyhow::bail!("invalid SV type: {:?}", arr[0]),
        };
        let chrom = arr[1].trim_start_matches("chr").to_string();
        let start = arr[2]
            .parse::<u32>()
            .map_err(|e| anyhow::anyhow!("could not parse start pos: {}", e))?;
        let stop = arr[3]
            .parse::<u32>()
            .map_err(|e| anyhow::anyhow!("could not parse stop pos: {}", e))?;

        Ok(StructuralVariant {
            chrom,
            start,
            stop,
            svtype,
            ambiguous_range: None,
        })
    }
}

impl Display for StructuralVariant {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}:{}:{}:{}",
            self.svtype, self.chrom, self.start, self.stop
        )
    }
}

impl From<StructuralVariant> for Interval {
    fn from(val: StructuralVariant) -> Self {
        Interval::new(
            val.chrom,
            val.start.saturating_sub(1).into()..val.stop.into(),
        )
    }
}

#[cfg(test)]
mod test {
    #[rstest::rstest]
    #[case("DEL:1:1000:2000")]
    #[case("DEL:chr1:1000:2000")]
    #[case("DUP:1:1000:2000")]
    #[case("DUP:chr1:1000:2000")]
    fn sv_from_str(#[case] val: &str) -> Result<(), anyhow::Error> {
        mehari::common::set_snapshot_suffix!("{}", val.replace(":", "-"));

        let sv = val.parse::<super::StructuralVariant>()?;
        insta::assert_yaml_snapshot!(sv);

        Ok(())
    }
}
