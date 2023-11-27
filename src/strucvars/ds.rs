//! Shared data structures for `strucvars`.

use std::{fmt::Display, str::FromStr};

use bio::bio_types::genome::Interval;

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

/// Mapping of gene identifiers.
#[derive(Debug, Clone, Default, PartialEq, serde::Deserialize, serde::Serialize)]
pub struct GeneIdInfo {
    /// HGNC identifier.
    pub hgnc_id: String,
    /// Official HGNC gene symbol.
    #[serde(alias = "gene_symbol")]
    pub hgnc_symbol: Option<String>,
    /// ENSEMBL gene identifier.
    pub ensembl_gene_id: Option<String>,
    /// NCBI gene identifier.
    #[serde(alias = "entrez_id")]
    pub ncbi_gene_id: Option<String>,
}

#[cfg(test)]
mod test {
    #[rstest::rstest]
    #[case("DEL:1:1000:2000")]
    #[case("DEL:chr1:1000:2000")]
    #[case("DUP:1:1000:2000")]
    #[case("DUP:chr1:1000:2000")]
    fn sv_from_str(#[case] val: &str) -> Result<(), anyhow::Error> {
        mehari::common::set_snapshot_suffix!("{}", val.replace(':', "-"));

        let sv = val.parse::<super::StructuralVariant>()?;
        insta::assert_yaml_snapshot!(sv);

        Ok(())
    }
}

pub mod intervals {
    //! Interval operations for rust-bio Intervals.

    use bio::bio_types::genome::AbstractInterval as _;

    /// The type to work with.
    pub type Interval = bio::bio_types::genome::Interval;

    /// Returns whether two intervals overlap.
    ///
    /// # Arguments
    ///
    /// * `lhs` - Left-hand side interval.
    /// * `rhs` - Right-hand side interval.
    ///
    /// # Returns
    ///
    /// `true` if the intervals overlap, `false` otherwise.
    pub fn do_overlap(lhs: &Interval, rhs: &Interval) -> bool {
        let lhs_contig = lhs.contig().strip_prefix("chr").unwrap_or(lhs.contig());
        let rhs_contig: &str = rhs.contig().strip_prefix("chr").unwrap_or(rhs.contig());
        if lhs_contig != rhs_contig {
            false
        } else {
            let a = lhs.range();
            let b = rhs.range();
            let intersect = std::cmp::max(a.start, b.start)..std::cmp::min(a.end, b.end);
            intersect.start < intersect.end
        }
    }

    /// Returns whether `lhs` interval contains the `rhs` one.
    ///
    /// # Arguments
    ///
    /// * `lhs` - Left-hand side interval.
    /// * `rhs` - Right-hand side interval.
    ///
    /// # Returns
    ///
    /// `true` if the intervals overlap, `false` otherwise.
    pub fn contains(lhs: &Interval, rhs: &Interval) -> bool {
        let lhs_contig = lhs.contig().strip_prefix("chr").unwrap_or(lhs.contig());
        let rhs_contig: &str = rhs.contig().strip_prefix("chr").unwrap_or(rhs.contig());
        if lhs_contig != rhs_contig {
            false
        } else {
            let a = lhs.range();
            let b = rhs.range();
            a.start <= b.start && a.end >= b.end
        }
    }

    /// Helper function that converts an `TxExonsRecord` into an Interval.
    pub fn exon_to_interval(
        chrom: String,
        record: &hgvs::data::interface::TxExonsRecord,
    ) -> Interval {
        Interval::new(
            chrom,
            (record.alt_start_i as u64)..(record.alt_end_i as u64),
        )
    }

    #[cfg(test)]
    mod test {
        use super::*;

        #[test]
        fn interval_overlaps_true() {
            let a = Interval::new("chr1".into(), 10..20);
            let b = Interval::new("chr1".into(), 15..25);
            assert!(do_overlap(&a, &b));
        }

        #[test]
        fn interval_overlaps_false() {
            let a = Interval::new("chr1".into(), 10..15);
            let b = Interval::new("chr1".into(), 15..25);
            assert!(!do_overlap(&a, &b));
        }

        #[test]
        fn interval_overlaps_different_contig() {
            let a = Interval::new("chr1".into(), 10..20);
            let b = Interval::new("chr2".into(), 15..25);
            assert!(!do_overlap(&a, &b));
        }

        #[test]
        fn interval_contains_true() {
            let a = Interval::new("chr1".into(), 10..20);
            let b = Interval::new("chr1".into(), 15..18);
            assert!(contains(&a, &b));
        }

        #[test]
        fn interval_contains_false() {
            let a = Interval::new("chr1".into(), 10..20);
            let b = Interval::new("chr1".into(), 15..25);
            assert!(!contains(&a, &b));
        }

        #[test]
        fn interval_contains_different_contig() {
            let a = Interval::new("chr1".into(), 10..20);
            let b = Interval::new("chr2".into(), 15..20);
            assert!(!contains(&a, &b));
        }
    }
}
