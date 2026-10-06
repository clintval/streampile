use std::ops::Range;

use super::{Read, bases, entries, entry, name, read, unfiltered};
use crate::auxiliary;
use crate::{ArraySubtype, AuxElement, AuxValue, EntryKind};

/// The views of the pileup of the reads at each position.
#[derive(Debug, PartialEq, Eq)]
struct Views {
    bases: String,
    qualities: Vec<u8>,
    unfiltered_depth: usize,
    filtered_depth: usize,
}

fn columns(reads: &[Read], positions: Range<usize>, floor: u8) -> Vec<Views> {
    let mut builder = unfiltered(reads).min_base_quality(floor);
    positions
        .map(|position| {
            let pileup = builder.pileup("chr1", position).expect("a pileup");
            Views {
                bases: bases(&pileup),
                qualities: pileup.qualities().collect(),
                unfiltered_depth: pileup.unfiltered_depth(),
                filtered_depth: pileup.filtered_depth(),
            }
        })
        .collect()
}

fn qualities(views: &[Views]) -> Vec<Vec<u8>> {
    views.iter().map(|views| views.qualities.clone()).collect()
}

fn base_lists(views: &[Views]) -> Vec<String> {
    views.iter().map(|views| views.bases.clone()).collect()
}

#[test]
fn test_a_base_entry_holds_its_base_and_quality() {
    let reads = [read("r", 10, "4M", "ACGT").quals(&[10, 20, 30, 40])];
    let mut builder = unfiltered(&reads);
    let pileup = builder.pileup("chr1", 12).unwrap();
    let entry = pileup.get(0).unwrap();
    assert_eq!((entry.base(), entry.quality()), (Some(b'G'), Some(30)));
    assert_eq!(
        (entry.query_position(), entry.query_position_or_next()),
        (Some(2), Some(2))
    );
}

#[test]
fn test_an_insertion_entry_holds_no_base_or_quality() {
    let reads = [read("r", 10, "2M2I2M", "ACGTAC").quals(&[10, 20, 30, 40, 50, 60])];
    let mut builder = unfiltered(&reads);
    let pileup = builder.pileup("chr1", 11).unwrap();
    let entry = pileup.get(1).unwrap();
    assert!(entry.is_insertion());
    assert_eq!((entry.base(), entry.quality()), (None, None));
    assert_eq!(
        entry.inserted_bases().map(Iterator::collect::<Vec<u8>>),
        Some(b"GT".to_vec())
    );
    assert_eq!(
        entry.inserted_qualities().map(Iterator::collect::<Vec<u8>>),
        Some(vec![30, 40])
    );
    assert_eq!(
        (entry.insertion_offset(), entry.insertion_length()),
        (Some(2), 2)
    );
}

#[test]
fn test_a_deletion_entry_holds_the_quality_of_the_next_base() {
    let reads = [read("r", 10, "1M1D3M", "ACGT").quals(&[10, 20, 30, 40])];
    let mut builder = unfiltered(&reads);
    let pileup = builder.pileup("chr1", 11).unwrap();
    let entry = pileup.get(0).unwrap();
    assert!(entry.is_deletion());
    assert_eq!((entry.base(), entry.quality()), (None, Some(20)));
    assert!(entry.inserted_bases().is_none() && entry.inserted_qualities().is_none());
}

#[test]
fn test_a_base_entry_of_an_n_is_a_no_call() {
    let reads = [read("r", 10, "4M", "AnGT"), read("s", 10, "1M2N1M", "AC")];
    let mut builder = unfiltered(&reads);
    let mut no_calls = Vec::new();
    for position in 10..14 {
        let pileup = builder.pileup("chr1", position).unwrap();
        no_calls.push(
            pileup
                .iter()
                .map(|entry| entry.is_no_call())
                .collect::<Vec<_>>(),
        );
    }
    assert_eq!(
        no_calls,
        [
            vec![false, false],
            vec![true, false],
            vec![false, false],
            vec![false, false]
        ]
    );
    let mut builder = unfiltered(&reads);
    let pileup = builder.pileup("chr1", 11).unwrap();
    let skip = pileup.get(1).unwrap();
    assert!(skip.is_skip() && !skip.is_no_call());
}

#[test]
fn test_an_entry_is_one_kind() {
    let reads = [read("r", 10, "1M1D1M1N1M1I1M", "ACGTA")];
    let mut builder = unfiltered(&reads).min_base_quality(0);
    let mut kinds = Vec::new();
    for position in 10..16 {
        let pileup = builder.pileup("chr1", position).unwrap();
        for entry in pileup.iter() {
            kinds.push((
                entry.kind(),
                entry.is_deletion(),
                entry.is_insertion(),
                entry.is_skip(),
            ));
        }
    }
    assert_eq!(
        kinds,
        [
            (EntryKind::Base, false, false, false),
            (EntryKind::Deletion, true, false, false),
            (EntryKind::Base, false, false, false),
            (EntryKind::Skip, false, false, true),
            (EntryKind::Base, false, false, false),
            (EntryKind::Insertion, false, true, false),
            (EntryKind::Base, false, false, false),
        ]
    );
}

#[test]
fn test_depths_count_bases_deletions_and_skips_but_not_insertions() {
    let base = read("base", 5, "10M", "AAAAAAAAAA");
    let deletion = read("deletion", 5, "3M3D4M", "AAATTTT");
    let skip = read("skip", 5, "3M3N4M", "AAATTTT");
    let insertion = read("insertion", 11, "1I6M", "AAATTTT");
    let cases: Vec<(Vec<Read>, (usize, usize))> = vec![
        (vec![], (0, 0)),
        (vec![base.clone()], (1, 1)),
        (vec![deletion.clone()], (1, 1)),
        (vec![skip.clone()], (1, 0)),
        (vec![insertion.clone()], (0, 0)),
        (vec![base, deletion, skip, insertion], (3, 2)),
    ];
    for (reads, depths) in cases {
        let mut builder = unfiltered(&reads);
        let pileup = builder.pileup("chr1", 10).unwrap();
        assert_eq!(
            (pileup.unfiltered_depth(), pileup.filtered_depth()),
            depths,
            "{reads:?}"
        );
    }
}

#[test]
fn test_qualities_of_one_read() {
    let reads = [read("r", 1, "4M", "ACGT").quals(&[20, 21, 22, 23])];
    assert_eq!(
        qualities(&columns(&reads, 1..5, 13)),
        [vec![20], vec![21], vec![22], vec![23]]
    );
    assert_eq!(
        qualities(&columns(&reads, 1..5, 22)),
        [vec![], vec![], vec![22], vec![23]]
    );
}

#[test]
fn test_qualities_leave_out_deletions_and_insertions() {
    let deleted = columns(
        &[read("r", 1, "2M1D1M", "ACG").quals(&[20, 21, 22])],
        1..5,
        13,
    );
    assert_eq!(qualities(&deleted), [vec![20], vec![21], vec![], vec![22]]);
    assert_eq!(deleted[2].filtered_depth, 1);
    let reads = [read("r", 1, "1M3I1M", "AGGGT").quals(&[31, 32, 33, 34, 35])];
    assert_eq!(
        qualities(&columns(&reads, 1..4, 13)),
        [vec![31], vec![35], vec![]]
    );
    let mut builder = unfiltered(&reads);
    let pileup = builder.pileup("chr1", 1).unwrap();
    assert_eq!(
        entries(&pileup)[1],
        entry("r", "insertion", None, None, Some("GGG"))
    );
}

#[test]
fn test_qualities_of_overlapping_reads() {
    let reads = [
        read("one", 1, "4M", "ACGT").quals(&[20, 21, 22, 23]),
        read("two", 3, "4M", "TGCA").quals(&[30, 31, 32, 33]),
    ];
    assert_eq!(
        qualities(&columns(&reads, 3..5, 13)),
        [vec![22, 30], vec![23, 31]]
    );
}

#[test]
fn test_bases_of_one_read() {
    let reads = [read("r", 1, "4M", "ACGT").quals(&[12, 13, 25, 30])];
    assert_eq!(
        base_lists(&columns(&reads, 0..6, 13)),
        ["", "", "C", "G", "T", ""]
    );
    assert_eq!(columns(&reads, 1..2, 0)[0].bases, "A");
}

#[test]
fn test_bases_leave_out_deletions_and_insertions() {
    let deleted = columns(&[read("r", 1, "2M1D1M", "ACG")], 1..5, 13);
    assert_eq!(base_lists(&deleted), ["A", "C", "", "G"]);
    let inserted = columns(&[read("r", 1, "1M3I1M", "AGGGT")], 1..4, 13);
    assert_eq!(base_lists(&inserted), ["A", "T", ""]);
}

#[test]
fn test_bases_of_overlapping_reads() {
    let reads = [read("one", 1, "4M", "ACGT"), read("two", 3, "4M", "TGCA")];
    assert_eq!(base_lists(&columns(&reads, 3..5, 13)), ["GT", "TG"]);
}

#[test]
fn test_views_of_reads_with_no_stored_bases_or_no_cigar_are_empty() {
    let no_bases = columns(&[read("r", 0, "4M", "*")], 0..1, 0);
    assert_eq!(
        no_bases[0],
        Views {
            bases: String::new(),
            qualities: vec![],
            unfiltered_depth: 1,
            filtered_depth: 0
        }
    );
    let mut builder = unfiltered(&[read("r", 0, "*", "*")]);
    assert!(builder.pileup("chr1", 0).unwrap().is_empty());
}

#[test]
fn test_without_overlaps_keeps_the_first_read_of_each_template() {
    let reads = [
        read("q3", 50, "50M", &"A".repeat(50)).flag(99),
        read("q1", 100, "50M", &"C".repeat(50)).flag(99),
        read("q2", 100, "50M", &"G".repeat(50)).flag(147),
        read("q3", 100, "50M", &"T".repeat(50)).flag(147),
        read("q1", 110, "50M", &"C".repeat(50)).flag(147),
        read("q2", 110, "50M", &"G".repeat(50)).flag(99),
    ];
    let mut builder = super::builder(&reads);
    assert_eq!(builder.pileup("chr1", 125).unwrap().unfiltered_depth(), 5);
    let mut builder = super::builder(&reads).without_overlaps(true);
    let kept = builder.pileup("chr1", 125).unwrap();
    assert_eq!(kept.unfiltered_depth(), 3);
    let kept_reads: Vec<(String, u16)> = kept
        .iter()
        .map(|entry| (name(entry.record()), entry.flags().bits()))
        .collect();
    assert_eq!(
        kept_reads,
        [
            ("q1".to_owned(), 99),
            ("q2".to_owned(), 147),
            ("q3".to_owned(), 147)
        ]
    );
    assert_eq!(bases(&kept), "CGT");
}

#[test]
fn test_without_overlaps_keeps_every_entry_of_the_kept_read() {
    let reads = [
        read("pair", 10, "3M2I3M", "ACGTTACG").flag(99),
        read("pair", 10, "6M", "ACGACG").flag(147),
        read("other", 12, "4M", "GACG"),
    ];
    let mut builder = unfiltered(&reads).min_base_quality(30);
    assert_eq!(builder.pileup("chr1", 12).unwrap().len(), 4);
    let mut builder = unfiltered(&reads)
        .min_base_quality(30)
        .without_overlaps(true);
    let kept = builder.pileup("chr1", 12).unwrap();
    assert_eq!(
        entries(&kept),
        [
            entry("pair", "base", Some(2), Some(2), None),
            entry("pair", "insertion", None, None, Some("TT")),
            entry("other", "base", Some(0), Some(0), None),
        ]
    );
    let flags: Vec<u16> = kept.iter().map(|entry| entry.flags().bits()).collect();
    assert_eq!(flags, [99, 99, 0]);
    assert_eq!(
        (kept.reference_sequence_name().to_string(), kept.position()),
        ("chr1".into(), 12)
    );
    assert_eq!(kept.min_base_quality(), 30);
}

#[test]
fn test_without_overlaps_keeps_a_mate_whose_base_another_mate_skips() {
    let mut bases_of_mate = "G".repeat(10);
    bases_of_mate.push('C');
    bases_of_mate.push_str(&"G".repeat(29));
    let reads = [
        read("t", 100, "20M300N20M", &"A".repeat(40)).flag(99),
        read("t", 330, "40M", &bases_of_mate).flag(147),
    ];
    let mut builder = unfiltered(&reads);
    let kinds: Vec<EntryKind> = builder
        .pileup("chr1", 340)
        .unwrap()
        .iter()
        .map(|entry| entry.kind())
        .collect();
    assert_eq!(kinds, [EntryKind::Skip, EntryKind::Base]);
    let mut builder = unfiltered(&reads).without_overlaps(true);
    let kept = builder.pileup("chr1", 340).unwrap();
    assert_eq!((kept.filtered_depth(), bases(&kept)), (1, "C".to_owned()));
    assert_eq!(kept.get(0).unwrap().flags().bits(), 147);
}

#[test]
fn test_without_overlaps_keeps_a_mate_at_the_floor_over_one_under_it() {
    let mate = read("t", 105, "10M", "GGCGGGGGGG").flag(147);
    let reads = [
        read("t", 100, "10M", &"A".repeat(10))
            .flag(99)
            .quals(&[2; 10]),
        mate.clone(),
    ];
    let mut builder = unfiltered(&reads).without_overlaps(true);
    let kept = builder.pileup("chr1", 107).unwrap();
    assert_eq!((kept.filtered_depth(), bases(&kept)), (1, "C".to_owned()));
    let reads = [reads[0].clone(), mate.quals(&[2; 10])];
    let mut builder = unfiltered(&reads).without_overlaps(true);
    let kept = builder.pileup("chr1", 107).unwrap();
    let flags: Vec<u16> = kept.iter().map(|entry| entry.flags().bits()).collect();
    assert_eq!((flags, kept.filtered_depth()), (vec![99], 0));
}

#[test]
fn test_without_overlaps_forgets_a_template_once_its_reads_are_passed() {
    let reads = [
        read("t", 10, "4M", "ACGT").flag(99),
        read("t", 12, "4M", "GTAC").flag(147),
        read("t", 100, "4M", "ACGT").flag(99),
    ];
    let mut builder = unfiltered(&reads).without_overlaps(true);
    let depths: Vec<usize> = [12, 15, 100]
        .into_iter()
        .map(|position| builder.pileup("chr1", position).unwrap().unfiltered_depth())
        .collect();
    assert_eq!(depths, [1, 1, 1]);
}

#[test]
fn test_pileups_leave_out_reads_on_other_contigs() {
    let reads = [
        read("a", 10, "4M", "ACGT"),
        read("b", 10, "4M", "TTTT").contig("chr2"),
    ];
    let mut builder = unfiltered(&reads);
    for (contig, read_name, base) in [("chr1", "a", "C"), ("chr2", "b", "T")] {
        let pileup = builder.pileup(contig, 11).unwrap();
        assert_eq!(super::names(&pileup), [read_name]);
        assert_eq!(bases(&pileup), base);
    }
}

#[test]
fn test_entries_read_aux_fields_at_their_query_offset() {
    let reads = [read("r", 10, "2M1D2M", "ACGT")
        .tag("ad:B:s,5,-6,7,8")
        .tag("ac:Z:ACGT")
        .tag("XF:f:1.5")
        .tag("XA:A:x")
        .tag("XI:i:-3")];
    for indexed in [false, true] {
        let mut builder = unfiltered(&reads);
        if indexed {
            builder = builder.index_aux_tags([*b"ad", *b"ac", *b"zz"]);
        }
        let pileup = builder.pileup("chr1", 11).unwrap();
        let entry = pileup.get(0).unwrap();
        assert_eq!(entry.aux_at(*b"ad").unwrap(), Some(AuxElement::Integer(-6)));
        assert_eq!(entry.aux_at(*b"ac").unwrap(), Some(AuxElement::Byte(b'C')));
        assert_eq!(entry.aux_at(*b"zz").unwrap(), None);
        assert_eq!(entry.aux(*b"XF").unwrap(), Some(AuxValue::Float(1.5)));
        assert_eq!(entry.aux(*b"XA").unwrap(), Some(AuxValue::Character(b'x')));
        assert_eq!(entry.aux(*b"XI").unwrap(), Some(AuxValue::Integer(-3)));
        let Some(AuxValue::Array(array)) = entry.aux(*b"ad").unwrap() else {
            panic!("ad is an array");
        };
        assert_eq!(
            (array.subtype(), array.len(), array.as_bytes().len()),
            (ArraySubtype::Int16, 4, 8)
        );
        let values: Vec<AuxElement> = array.iter().collect();
        assert_eq!(values, [5, -6, 7, 8].map(AuxElement::Integer));
        assert_eq!(array.get(4), None);
        let pileup = builder.pileup("chr1", 12).unwrap();
        let deletion = pileup.get(0).unwrap();
        assert!(deletion.is_deletion());
        assert_eq!(deletion.aux_at(*b"ad").unwrap(), None);
        let pileup = builder.pileup("chr1", 13).unwrap();
        assert_eq!(
            pileup.get(0).unwrap().aux_at(*b"ad").unwrap(),
            Some(AuxElement::Integer(7))
        );
    }
}

#[test]
fn test_malformed_aux_fields_are_errors() {
    for data in [
        &b"XAZab"[..],
        b"XA",
        b"XAB",
        b"XABs\xff\xff\xff\x7f",
        b"XABq\x01\x00\x00\x00\x00",
        b"XAi\x01\x00",
        b"XAQ\x01",
    ] {
        assert!(
            auxiliary::find(data, *b"XA").is_err(),
            "{}",
            data.escape_ascii()
        );
    }
    assert_eq!(auxiliary::find(b"", *b"XA").unwrap(), None);
    let truncated_after = b"XAc\x01XBi\x01";
    assert_eq!(
        auxiliary::find(truncated_after, *b"XA").unwrap(),
        Some(AuxValue::Integer(1))
    );
    assert!(auxiliary::find(truncated_after, *b"XB").is_err());
    let data = b"XAc\xffXBC\xffXCs\xfe\xffXDS\xfe\xffXEI\x01\x00\x00\x80XFHff\x00";
    let found: Vec<Option<AuxValue<'_>>> = [*b"XA", *b"XB", *b"XC", *b"XD", *b"XE", *b"XF"]
        .map(|tag| auxiliary::find(data, tag).unwrap())
        .into();
    assert_eq!(
        found,
        [
            Some(AuxValue::Integer(-1)),
            Some(AuxValue::Integer(255)),
            Some(AuxValue::Integer(-2)),
            Some(AuxValue::Integer(65_534)),
            Some(AuxValue::Integer(2_147_483_649)),
            Some(AuxValue::Hex(b"ff")),
        ]
    );
}

#[test]
fn test_every_array_subtype_is_read() {
    let arrays: [(&[u8], AuxElement); 7] = [
        (b"Bc\x01\x00\x00\x00\xff", AuxElement::Integer(-1)),
        (b"BC\x01\x00\x00\x00\xff", AuxElement::Integer(255)),
        (b"Bs\x01\x00\x00\x00\xff\xff", AuxElement::Integer(-1)),
        (b"BS\x01\x00\x00\x00\xff\xff", AuxElement::Integer(65_535)),
        (
            b"Bi\x01\x00\x00\x00\xff\xff\xff\xff",
            AuxElement::Integer(-1),
        ),
        (
            b"BI\x01\x00\x00\x00\xff\xff\xff\xff",
            AuxElement::Integer(4_294_967_295),
        ),
        (
            b"Bf\x01\x00\x00\x00\x00\x00\xc0\x3f",
            AuxElement::Float(1.5),
        ),
    ];
    for (value, element) in arrays {
        let data = [&b"XA"[..], value].concat();
        let Some(AuxValue::Array(array)) = auxiliary::find(&data, *b"XA").unwrap() else {
            panic!("{} is an array", value.escape_ascii());
        };
        assert_eq!(
            (array.len(), array.is_empty(), array.get(0)),
            (1, false, Some(element))
        );
        assert_eq!(AuxValue::Array(array).get(0), Some(element));
    }
    assert_eq!(AuxValue::Integer(1).get(0), None);
}
