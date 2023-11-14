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
    if lhs.contig() != rhs.contig() {
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
    if lhs.contig() != rhs.contig() {
        false
    } else {
        let a = lhs.range();
        let b = rhs.range();
        a.start <= b.start && a.end >= b.end
    }
}

/// Helper function that converts an `TxExonsRecord` into an Interval.
pub fn exon_to_interval(chrom: String, record: &hgvs::data::interface::TxExonsRecord) -> Interval {
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
