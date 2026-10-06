use super::{HEADER, read, records};
use crate::footprint::{Footprint, Insertion, Located};

const SKIPPED: i64 = -1;
const DELETED_AT_END: i64 = -2;

/// The offset Python's `Footprint` stores for a deletion followed by the base at `offset`.
fn next_after(offset: i64) -> i64 {
    -offset - 3
}

fn located(offset: i64) -> Located {
    match offset {
        SKIPPED => Located::Skip,
        DELETED_AT_END => Located::Deletion(None),
        offset if offset >= 0 => Located::Base(offset as u32),
        offset => Located::Deletion(Some((-offset - 3) as u32)),
    }
}

fn footprint(cigar: &str, bases: &str) -> (Footprint, bool) {
    let (_, records) = records(HEADER, &[read("r", 10, cigar, bases)]);
    let record = &records[0];
    let mut footprint = Footprint::default();
    let placed = footprint
        .fill(10, record.cigar().as_bytes(), record.sequence().len())
        .expect("a valid CIGAR");
    (footprint, placed)
}

/// A CIGAR, its bases, the offsets Python's `Footprint` stores, its insertions, and its end.
type Case = (
    &'static str,
    &'static str,
    Vec<i64>,
    Vec<(i64, u32, u32)>,
    i64,
);

#[test]
#[allow(clippy::too_many_lines)]
fn test_footprint_walks_the_cigar_once() {
    let cases: Vec<Case> = vec![
        ("4M", "ACGT", vec![0, 1, 2, 3], vec![], 14),
        ("2S4M", "TTACGT", vec![2, 3, 4, 5], vec![], 14),
        ("3H4M2H", "ACGT", vec![0, 1, 2, 3], vec![], 14),
        (
            "2M2D2M",
            "ACGT",
            vec![0, 1, next_after(2), next_after(2), 2, 3],
            vec![],
            16,
        ),
        (
            "2M2N2M",
            "ACGT",
            vec![0, 1, SKIPPED, SKIPPED, 2, 3],
            vec![],
            16,
        ),
        ("2M1I2M", "ACGTA", vec![0, 1, 3, 4], vec![(11, 2, 1)], 14),
        ("1I3M", "ACGT", vec![1, 2, 3], vec![(9, 0, 1)], 13),
        ("2S1I3M", "TTACGT", vec![3, 4, 5], vec![(9, 2, 1)], 13),
        ("3M1I", "ACGT", vec![0, 1, 2], vec![(12, 3, 1)], 13),
        (
            "2M1D1I2M",
            "ACGTA",
            vec![0, 1, next_after(2), 3, 4],
            vec![(12, 2, 1)],
            15,
        ),
        ("2M1I1I2M", "ACGTAC", vec![0, 1, 4, 5], vec![(11, 2, 2)], 14),
        (
            "3M2D",
            "ACG",
            vec![0, 1, 2, DELETED_AT_END, DELETED_AT_END],
            vec![],
            15,
        ),
        (
            "3M2D1S",
            "ACGT",
            vec![0, 1, 2, next_after(3), next_after(3)],
            vec![],
            15,
        ),
        (
            "3M2D1M",
            "*",
            vec![0, 1, 2, next_after(3), next_after(3), 3],
            vec![],
            16,
        ),
        (
            "3M2D",
            "*",
            vec![0, 1, 2, DELETED_AT_END, DELETED_AT_END],
            vec![],
            15,
        ),
        ("1=1X2M", "ACGT", vec![0, 1, 2, 3], vec![], 14),
        ("1D3M", "ACG", vec![next_after(0), 0, 1, 2], vec![], 14),
        (
            "3D4M",
            "ACGT",
            vec![next_after(0), next_after(0), next_after(0), 0, 1, 2, 3],
            vec![],
            17,
        ),
        (
            "1D1I3M",
            "ACGT",
            vec![next_after(0), 1, 2, 3],
            vec![(10, 0, 1)],
            14,
        ),
        ("4H1I3M", "ACGT", vec![1, 2, 3], vec![(9, 0, 1)], 13),
        ("4M4D", "ACGT", vec![0, 1, 2, 3, -2, -2, -2, -2], vec![], 18),
        (
            "2M1I1P1I2M",
            "ACGTAC",
            vec![0, 1, 4, 5],
            vec![(11, 2, 2)],
            14,
        ),
    ];
    for (cigar, bases, offsets, insertions, end) in cases {
        let (mut footprint, placed) = footprint(cigar, bases);
        assert!(placed, "{cigar}");
        assert_eq!((footprint.start, footprint.end), (10, end), "{cigar}");
        let walked: Vec<Located> = (footprint.start..footprint.end)
            .filter_map(|pos| footprint.locate(pos))
            .collect();
        let expected: Vec<Located> = offsets.into_iter().map(located).collect();
        assert_eq!(walked, expected, "{cigar}");
        let (mut footprint, _) = self::footprint(cigar, bases);
        let found: Vec<Insertion> = (0..30)
            .filter_map(|pos| footprint.insertion_at(pos))
            .collect();
        let expected: Vec<Insertion> = insertions
            .into_iter()
            .map(|(anchor, offset, length)| Insertion {
                anchor,
                offset,
                length,
            })
            .collect();
        assert_eq!(found, expected, "{cigar}");
    }
}

#[test]
fn test_footprint_locates_nothing_outside_the_read() {
    let (mut footprint, _) = footprint("2S4M", "TTACGT");
    assert_eq!(footprint.locate(9), None);
    assert_eq!(footprint.locate(10), Some(Located::Base(2)));
    assert_eq!(footprint.locate(14), None);
}

#[test]
fn test_a_long_reference_skip_is_one_block() {
    let (mut footprint, _) = footprint("2M100000N2M", "ACGT");
    assert_eq!(footprint.end, 10 + 100_004);
    assert_eq!(footprint.locate(50_000), Some(Located::Skip));
    assert_eq!(footprint.locate(100_013), Some(Located::Base(3)));
}

#[test]
fn test_footprint_is_placed_only_with_a_reference_consuming_operator() {
    for (cigar, bases, placed) in [
        ("4M", "ACGT", true),
        ("1D3M", "ACG", true),
        ("2S2N2M", "ACGT", true),
        ("1=1X2S", "ACGT", true),
        ("4D", "*", true),
        ("4S", "ACGT", false),
        ("4I", "ACGT", false),
        ("2S2I", "ACGT", false),
        ("4H4S", "ACGT", false),
        ("1S2I1S", "ACGT", false),
    ] {
        assert_eq!(footprint(cigar, bases).1, placed, "{cigar}");
    }
}

#[test]
fn test_footprint_refuses_a_cigar_with_a_partial_operator() {
    let mut footprint = Footprint::default();
    assert!(footprint.fill(0, &[0x40, 0, 0], 4).is_err());
    assert!(footprint.fill(0, &[0x49, 0, 0, 0], 4).is_err());
}

#[test]
fn test_footprint_refuses_a_cigar_whose_query_length_differs_from_the_sequence() {
    let four_matches = (4_u32 << 4).to_le_bytes();
    let mut footprint = Footprint::default();
    for stored in [3, 5] {
        let refused = footprint.fill(10, &four_matches, stored).unwrap_err();
        assert_eq!(
            refused.to_string(),
            "CIGAR and query sequence lengths differ"
        );
    }
    assert!(footprint.fill(10, &four_matches, 4).unwrap());
    assert!(footprint.fill(10, &four_matches, 0).unwrap());
    assert_eq!(footprint.query_length, 4);
}
