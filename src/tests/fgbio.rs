//! Tests that carry over the intent of fgbio's pileup tests, with fgbio's 1-based positions
//! moved to 0-based ones.

use std::collections::BTreeMap;

use super::{Read, bases, builder, name, names, read};
use crate::builder::mate_span;
use crate::{AuxValue, EntryKind, Error};

const READ_LENGTH: usize = 50;

fn repeat(base: char) -> String {
    base.to_string().repeat(READ_LENGTH)
}

/// A pair of 50-base reads at 0-based starts with mate fields as htsjdk sets them, sorted.
fn pair(name: &str, start1: usize, start2: usize, reverse1: bool, reverse2: bool) -> Vec<Read> {
    let end = |start: usize| start + READ_LENGTH - 1;
    let five_prime = |start: usize, reverse: bool| if reverse { end(start) } else { start } as i64;
    let (first, second) = (five_prime(start1, reverse1), five_prime(start2, reverse2));
    let insert = (second - first + if second >= first { 1 } else { -1 }) as i32;
    let strand =
        |reverse: bool, mate_reverse: bool| u16::from(reverse) * 16 + u16::from(mate_reverse) * 32;
    let cigar = format!("{READ_LENGTH}M");
    let r1 = read(name, start1, &cigar, &repeat('A'))
        .flag(1 | 64 | strand(reverse1, reverse2))
        .mate("chr1", start2, insert)
        .tag(&format!("MC:Z:{cigar}"));
    let r2 = read(name, start2, &cigar, &repeat('C'))
        .flag(1 | 128 | strand(reverse2, reverse1))
        .mate("chr1", start1, -insert)
        .tag(&format!("MC:Z:{cigar}"));
    if start2 < start1 {
        vec![r2, r1]
    } else {
        vec![r1, r2]
    }
}

fn counts(bases: &str) -> BTreeMap<char, usize> {
    let mut counts = BTreeMap::new();
    for base in bases.chars() {
        *counts.entry(base).or_default() += 1;
    }
    counts
}

fn sorted(bases: &str) -> String {
    let mut bases: Vec<char> = bases.chars().collect();
    bases.sort_unstable();
    bases.into_iter().collect()
}

#[test]
fn test_builder_leaves_out_reads_on_a_previous_contig_or_ending_just_before() {
    let mut reads = Vec::new();
    for _ in 0..5 {
        reads.push(read("prev-contig", 54, "50M", &repeat('A')));
    }
    for _ in 0..5 {
        reads.push(read("abutting", 4, "50M", &repeat('A')).contig("chr2"));
    }
    for (base, count) in [('A', 5), ('C', 4), ('G', 3), ('T', 2), ('N', 1)] {
        for _ in 0..count {
            reads.push(read("here", 54, "50M", &repeat(base)).contig("chr2"));
        }
    }
    let mut builder = builder(&reads);
    let pileup = builder.pileup("chr2", 54).unwrap();
    assert_eq!(pileup.unfiltered_depth(), 15);
    assert!(names(&pileup).iter().all(|read_name| read_name == "here"));
    let expected = BTreeMap::from([('A', 5), ('C', 4), ('G', 3), ('N', 1), ('T', 2)]);
    assert_eq!(counts(&bases(&pileup)), expected);
}

#[test]
fn test_builder_piles_up_every_edge_case_of_indels() {
    let reads = [
        read("q1", 100, "10M2D40M", &repeat('A')),
        read("q2", 100, "10M2I38M", &repeat('C')),
        read("q3", 100, "31M9I10M", &repeat('G')),
        read("q4", 100, "30M9D20M", &repeat('T')),
        read("q5", 140, "10I40M", &repeat('N')),
        read("q6", 200, "47M3S", &repeat('N')),
    ];
    let mut builder = builder(&reads);
    let named = |kind: EntryKind, pileup: &crate::Pileup<'_>| -> Vec<String> {
        pileup
            .iter()
            .filter(|entry| entry.kind() == kind)
            .map(|entry| name(entry.record()))
            .collect()
    };

    assert_eq!(bases(&builder.pileup("chr1", 104).unwrap()).len(), 4);

    let before = builder.pileup("chr1", 109).unwrap();
    assert_eq!((before.unfiltered_depth(), before.len()), (4, 5));
    assert_eq!(sorted(&bases(&before)), "ACGT");
    assert_eq!(named(EntryKind::Insertion, &before), ["q2"]);

    let deleted = builder.pileup("chr1", 110).unwrap();
    assert_eq!((deleted.unfiltered_depth(), deleted.len()), (4, 4));
    assert_eq!(sorted(&bases(&deleted)), "CGT");
    assert_eq!(named(EntryKind::Deletion, &deleted), ["q1"]);

    let bigger = builder.pileup("chr1", 130).unwrap();
    assert_eq!((bigger.unfiltered_depth(), bigger.len()), (4, 5));
    assert_eq!(sorted(&bases(&bigger)), "ACG");
    assert_eq!(named(EntryKind::Insertion, &bigger), ["q3"]);
    assert_eq!(named(EntryKind::Deletion, &bigger), ["q4"]);

    let leading = builder.pileup("chr1", 139).unwrap();
    assert_eq!((leading.unfiltered_depth(), leading.len()), (4, 5));
    assert_eq!(sorted(&bases(&leading)), "ACGT");
    assert_eq!(named(EntryKind::Insertion, &leading), ["q5"]);
    assert_eq!(leading.get(4).unwrap().insertion_offset(), Some(0));

    let clipped = builder.pileup("chr1", 246).unwrap();
    assert_eq!(
        (clipped.unfiltered_depth(), clipped.len(), bases(&clipped)),
        (1, 1, "N".into())
    );

    let past = builder.pileup("chr1", 247).unwrap();
    assert_eq!((past.unfiltered_depth(), past.len()), (0, 0));
}

#[test]
fn test_builder_piles_up_only_reads_of_mapped_pairs_with_a_read_filter() {
    let mut half_mapped = pair("q2", 100, 100, false, false);
    half_mapped[1] = half_mapped[1].clone().flag(1 | 4 | 128);
    half_mapped[0] = half_mapped[0].clone().flag(1 | 8 | 64);
    let mut reads = vec![read("q1", 100, "50M", &repeat('A'))];
    reads.extend(half_mapped);
    reads.extend(pair("q3", 100, 299, false, true));
    let mapped_pair = |record: &noodles::bam::Record| {
        let flags = record.flags();
        flags.is_segmented() && !flags.is_unmapped() && !flags.is_mate_unmapped()
    };
    let mut builder = builder(&reads).read_filter(mapped_pair);
    let pileup = builder.pileup("chr1", 104).unwrap();
    assert_eq!(pileup.unfiltered_depth(), 1);
    let kept = pileup.get(0).unwrap();
    assert_eq!(
        (name(kept.record()), kept.flags().is_first_segment()),
        ("q3".into(), true)
    );
}

#[test]
fn test_builder_keeps_a_fragment_and_positions_outside_the_insert_of_an_fr_pair() {
    let mut fragment = builder(&[read("q1", 99, "50M", &repeat('A'))]);
    assert_eq!(fragment.pileup("chr1", 99).unwrap().unfiltered_depth(), 1);
    assert_eq!(fragment.pileup("chr1", 148).unwrap().unfiltered_depth(), 1);
    let pileup = fragment.pileup("chr1", 148).unwrap();
    assert_eq!(
        pileup.get(0).unwrap().template_end_distance().unwrap(),
        None
    );

    let mut builder = builder(&pair("q2", 100, 99, false, true));
    let mut depths = Vec::new();
    let mut outside = Vec::new();
    for position in [99, 100, 148, 149] {
        let pileup = builder.pileup("chr1", position).unwrap();
        depths.push(pileup.unfiltered_depth());
        outside.push(
            pileup
                .iter()
                .filter(|entry| entry.template_end_distance().unwrap().is_none())
                .count(),
        );
    }
    assert_eq!(depths, [1, 2, 2, 1]);
    assert_eq!(outside, [1, 0, 0, 1]);
}

#[test]
fn test_builder_keeps_every_position_of_a_pair_with_the_reverse_read_starting_later() {
    let mut builder = builder(&pair("q2", 100, 99, true, false));
    let mut depths = Vec::new();
    let mut distances = Vec::new();
    for position in [99, 100, 148, 149] {
        let pileup = builder.pileup("chr1", position).unwrap();
        depths.push(pileup.unfiltered_depth());
        distances.extend(
            pileup
                .iter()
                .map(|entry| entry.template_end_distance().unwrap()),
        );
    }
    assert_eq!(depths, [1, 2, 2, 1]);
    let expected = [50, 49, 1, 1, 49, 50].map(Some);
    assert_eq!(distances, expected);
}

#[test]
fn test_builder_composes_a_read_filter_with_an_entry_filter() {
    let reads: Vec<Read> = [("q1", 100), ("x2", 104), ("q3", 108), ("q4", 112)]
        .map(|(read_name, start)| read(read_name, start, "50M", &repeat('A')))
        .into();
    let mut builder = builder(&reads).read_filter(|record| !name(record).starts_with('x'));
    let pileup = builder.pileup("chr1", 114).unwrap();
    let kept: Vec<String> = pileup
        .iter()
        .filter(|entry| entry.query_position().is_some_and(|offset| offset > 5))
        .map(|entry| name(entry.record()))
        .collect();
    assert_eq!(kept, ["q1", "q3"]);
}

#[test]
fn test_entries_report_offsets_and_positions_in_read_order() {
    let mut reads = pair("q1", 100, 200, false, true);
    reads[0] = reads[0].clone().quals(&[35; READ_LENGTH]);
    reads[1] = reads[1].clone().quals(&[35; READ_LENGTH]);
    let mut builder = builder(&reads);
    let mut seen = Vec::new();
    for position in [104, 204] {
        let pileup = builder.pileup("chr1", position).unwrap();
        assert_eq!(pileup.unfiltered_depth(), 1);
        let entry = pileup.get(0).unwrap();
        seen.push((
            entry.base(),
            entry.sequenced_base(),
            entry.quality(),
            entry.query_position(),
            entry.five_prime_distance(),
        ));
    }
    assert_eq!(
        seen,
        [
            (Some(b'A'), Some(b'A'), Some(35), Some(4), Some(4)),
            (Some(b'C'), Some(b'G'), Some(35), Some(4), Some(45)),
        ]
    );
}

#[test]
fn test_a_template_end_needs_a_mapped_fr_mate_on_the_same_contig() {
    let mut unmapped_mate = pair("q2", 100, 100, false, false);
    unmapped_mate[0] = unmapped_mate[0].clone().flag(1 | 8 | 64);
    unmapped_mate[1] = unmapped_mate[1].clone().flag(1 | 4 | 128);
    let other_contig = read("q3", 100, "50M", &repeat('A'))
        .flag(1 | 32 | 64)
        .mate("chr2", 200, 0)
        .tag("MC:Z:50M");
    let cases = [
        vec![read("q1", 100, "50M", &repeat('A'))],
        unmapped_mate,
        vec![other_contig],
        pair("pp", 100, 200, false, false),
        pair("mm", 100, 200, true, true),
        pair("rf", 100, 200, true, false),
    ];
    for reads in cases {
        let mut builder = builder(&reads);
        for position in [100, 149, 200, 249] {
            let pileup = builder.pileup("chr1", position).unwrap();
            assert!(
                pileup
                    .iter()
                    .all(|entry| entry.template_end_distance().unwrap().is_none()),
                "{reads:?}"
            );
        }
    }
}

#[test]
fn test_a_template_end_spans_the_insert_of_an_fr_pair() {
    let short = |reads: Vec<Read>| -> Vec<Read> {
        reads
            .into_iter()
            .map(|read| Read {
                cigar: "10M".into(),
                bases: "A".repeat(10),
                quals: Some(vec![40; 10]),
                tags: vec!["MC:Z:10M".into()],
                ..read
            })
            .collect()
    };
    let mut reads = short(pair("q", 99, 190, false, true));
    reads[0].tlen = 101;
    reads[1].tlen = -101;
    let mut builder = builder(&reads);
    let first = builder
        .pileup("chr1", 99)
        .unwrap()
        .get(0)
        .unwrap()
        .template_end_distance()
        .unwrap();
    let last = builder
        .pileup("chr1", 199)
        .unwrap()
        .get(0)
        .unwrap()
        .template_end_distance()
        .unwrap();
    assert_eq!((first, last), (Some(100), Some(100)));
}

#[test]
fn test_a_template_end_is_measured_from_the_mate_cigar_alone() {
    let expected = [
        (false, Some(99)),
        (false, Some(89)),
        (false, Some(50)),
        (true, Some(50)),
        (true, Some(89)),
        (true, Some(99)),
    ];
    let with_mate_cigar = pair("q", 100, 150, false, true);
    let mut any_template_length = with_mate_cigar.clone();
    for (read, tlen) in any_template_length.iter_mut().zip([7, 0]) {
        read.tlen = tlen;
    }
    let without: Vec<Read> = with_mate_cigar
        .iter()
        .cloned()
        .map(|read| Read {
            tags: vec![],
            ..read
        })
        .collect();
    for reads in [with_mate_cigar, any_template_length, without] {
        let has_mate_cigar = !reads[0].tags.is_empty();
        let mut builder = builder(&reads);
        let mut distances = Vec::new();
        for position in [100, 110, 149, 150, 189, 199] {
            let pileup = builder.pileup("chr1", position).unwrap();
            for entry in pileup.iter() {
                match entry.template_end_distance() {
                    Ok(distance) => distances.push((entry.is_reverse(), distance)),
                    Err(Error::MissingMateCigar { name }) => {
                        assert!(!has_mate_cigar && !entry.is_reverse() && name == "q");
                    }
                    Err(error) => panic!("{error}"),
                }
            }
        }
        let expected = if has_mate_cigar {
            &expected[..]
        } else {
            &expected[3..]
        };
        assert_eq!(distances, expected, "{reads:?}");
    }
}

#[test]
fn test_a_template_end_from_an_invalid_mate_cigar_is_an_error_naming_the_read() {
    for value in ["10M5", "4Q", "M", "0M", "2S", "*", "4294967295M1M"] {
        let mut reads = pair("q", 100, 150, false, true);
        reads[0].tags = vec![format!("MC:Z:{value}")];
        let mut builder = builder(&reads);
        let pileup = builder.pileup("chr1", 120).unwrap();
        assert_eq!(pileup.len(), 1);
        let refused = pileup.get(0).unwrap().template_end_distance();
        assert!(
            matches!(
                &refused,
                Err(Error::InvalidMateCigar { name, value: found }) if name == "q" && found == value
            ),
            "{value}: {refused:?}"
        );
    }
    let mut reads = pair("q", 100, 150, false, true);
    reads[0].tags = vec!["MC:i:50".into()];
    let mut builder = builder(&reads);
    let refused = builder
        .pileup("chr1", 120)
        .unwrap()
        .get(0)
        .unwrap()
        .template_end_distance();
    assert!(matches!(refused, Err(Error::InvalidMateCigar { .. })));
}

#[test]
fn test_the_span_of_a_mate_cigar() {
    let span = |text: &str| mate_span(AuxValue::String(text.as_bytes()));
    assert_eq!(span("10M2I3D4N1=1X5S2H1P"), Some(19));
    assert_eq!(span("+3M"), Some(3));
    for invalid in [
        "",
        "M",
        "0M",
        "2S",
        "10",
        "10Q",
        "99999999999999999999M",
        "9223372036854775807M9223372036854775807M",
        "268435456M",
    ] {
        assert_eq!(span(invalid), None, "{invalid}");
    }
    assert_eq!(span("268435455M268435455N"), Some(536_870_910));
    assert_eq!(mate_span(AuxValue::Integer(10)), None);
}
