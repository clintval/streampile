//! Tests that carry over the intent of fgbio's pileup tests, with fgbio's 1-based positions
//! moved to 0-based ones. Pairs are built at fgbio's 1-based starts by the shared `SamBuilder`,
//! as fgbio builds them, and other reads at 0-based starts.

use std::collections::BTreeMap;

use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::record_buf::data::field::value::Array;

use super::{Read, bases, builder, name, names, read};
use crate::template::{Alignment, parse_cigar};
use crate::testing::{BuiltRecords, Frag, Pair, SamBuilder, Strand};
use crate::{EntryKind, Error, StreamingPileupBuilder};

const READ_LENGTH: usize = 50;

fn repeat(base: char) -> String {
    base.to_string().repeat(READ_LENGTH)
}

/// A builder of fgbio's 50-base reads.
fn reads() -> SamBuilder {
    SamBuilder::new().read_length(READ_LENGTH)
}

/// A pair of 50-base reads of `A`s and `C`s at 1-based starts on these strands.
fn pair(
    name: &str,
    start1: usize,
    start2: usize,
    strand1: Strand,
    strand2: Strand,
) -> Vec<RecordBuf> {
    reads().add_pair(Pair {
        name: Some(name.into()),
        bases1: Some(repeat('A')),
        bases2: Some(repeat('C')),
        strand1,
        strand2,
        ..Pair::at(start1, start2)
    })
}

/// A pileup builder over built records, changed or not, in coordinate order.
fn piled(
    records: impl IntoIterator<Item = RecordBuf>,
) -> StreamingPileupBuilder<'static, BuiltRecords> {
    let mut builder = reads();
    builder.extend(records);
    builder.to_pileup_builder()
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
    let mut builder = reads();
    builder.add_frag(Frag {
        name: Some("q1".into()),
        ..Frag::at(101)
    });
    builder.add_pair(Pair {
        name: Some("q2".into()),
        unmapped2: true,
        ..Pair::at(101, 101)
    });
    builder.add_pair(Pair {
        name: Some("q3".into()),
        ..Pair::at(101, 300)
    });
    let mapped_pair = |record: &noodles::bam::Record| {
        let flags = record.flags();
        flags.is_segmented() && !flags.is_unmapped() && !flags.is_mate_unmapped()
    };
    let mut builder = builder.to_pileup_builder().read_filter(mapped_pair);
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
    let mut fragment = reads();
    fragment.add_frag(Frag::at(100));
    let mut fragment = fragment.to_pileup_builder();
    assert_eq!(fragment.pileup("chr1", 99).unwrap().unfiltered_depth(), 1);
    assert_eq!(fragment.pileup("chr1", 148).unwrap().unfiltered_depth(), 1);
    let pileup = fragment.pileup("chr1", 148).unwrap();
    assert_eq!(
        pileup.get(0).unwrap().template_end_distance().unwrap(),
        None
    );

    let mut builder = piled(pair("q2", 101, 100, Strand::Plus, Strand::Minus));
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
    let mut builder = piled(pair("q2", 101, 100, Strand::Minus, Strand::Plus));
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
    let mut builder = reads().base_quality(35);
    builder.add_pair(Pair {
        name: Some("q1".into()),
        bases1: Some(repeat('A')),
        bases2: Some(repeat('C')),
        ..Pair::at(101, 201)
    });
    let mut builder = builder.to_pileup_builder();
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
    let fragment = reads().add_frag(Frag::at(101));
    let unmapped_mate = reads().add_pair(Pair {
        unmapped2: true,
        ..Pair::at(101, 101)
    });
    let other_contig = reads().add_pair(Pair {
        contig2: Some(1),
        ..Pair::at(101, 201)
    });
    let cases = [
        fragment,
        unmapped_mate,
        other_contig,
        pair("pp", 101, 201, Strand::Plus, Strand::Plus),
        pair("mm", 101, 201, Strand::Minus, Strand::Minus),
        pair("rf", 101, 201, Strand::Minus, Strand::Plus),
    ];
    for records in cases {
        let mut builder = piled(records.clone());
        for position in [100, 149, 200, 249] {
            let pileup = builder.pileup("chr1", position).unwrap();
            assert!(
                pileup
                    .iter()
                    .all(|entry| entry.template_end_distance().unwrap().is_none()),
                "{records:?}"
            );
        }
    }
}

#[test]
fn test_a_template_end_spans_the_insert_of_an_fr_pair() {
    let records = SamBuilder::new()
        .read_length(10)
        .add_pair(Pair::at(100, 191));
    assert_eq!(records[0].template_length(), 101);
    let mut builder = piled(records);
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
    let with_mate_cigar = reads().add_pair(Pair {
        name: Some("q".into()),
        ..Pair::at(101, 151)
    });
    let mut any_template_length = with_mate_cigar.clone();
    for (record, tlen) in any_template_length.iter_mut().zip([7, 0]) {
        *record.template_length_mut() = tlen;
    }
    let without: Vec<RecordBuf> = with_mate_cigar
        .iter()
        .cloned()
        .map(SamBuilder::without_mate_cigar)
        .collect();
    for records in [with_mate_cigar, any_template_length, without] {
        let has_mate_cigar = records[0].data().get(&Tag::MATE_CIGAR).is_some();
        let mut builder = piled(records.clone());
        let mut distances = Vec::new();
        for position in [100, 110, 149, 150, 189, 199] {
            let pileup = builder.pileup("chr1", position).unwrap();
            for entry in pileup.iter() {
                match entry.template_end_distance() {
                    Ok(distance) => distances.push((entry.is_reverse(), distance)),
                    Err(Error::MissingMateCigar { name }) => {
                        assert!(!has_mate_cigar && name == "q");
                    }
                    Err(error) => panic!("{error}"),
                }
            }
        }
        let expected = if has_mate_cigar { &expected[..] } else { &[] };
        assert_eq!(distances, expected, "{records:?}");
    }
}

#[test]
fn test_a_template_end_from_an_invalid_mate_cigar_is_an_error_naming_the_read() {
    let with_mate_cigar = |value: Value| {
        let mut records = pair("q", 101, 151, Strand::Plus, Strand::Minus);
        records[0].data_mut().insert(Tag::MATE_CIGAR, value);
        let mut builder = piled(records);
        let pileup = builder.pileup("chr1", 120).unwrap();
        assert_eq!(pileup.len(), 1);
        pileup.get(0).unwrap().template_end_distance()
    };
    for value in ["10M5", "4Q", "M", "0M", "2S", "*", "4294967295M1M"] {
        let refused = with_mate_cigar(Value::from(value));
        assert!(
            matches!(
                &refused,
                Err(Error::InvalidMateCigar { name, value: found }) if name == "q" && found == value
            ),
            "{value}: {refused:?}"
        );
    }
    for (value, shown) in [
        (Value::from(50_i32), "50"),
        (Value::from(4_000_000_000_u32), "4000000000"),
        (Value::Float(1.5), "1.5"),
        (Value::Character(b'M'), "M"),
        (Value::Array(Array::Int32(vec![1, -2])), "1,-2"),
    ] {
        let refused = with_mate_cigar(value).unwrap_err();
        assert_eq!(
            refused.to_string(),
            format!("read q has an invalid MC tag: {shown}")
        );
    }
}

/// A forward read and its reverse mate at fgbio's 1-based starts, as fgbio's `SamBuilder.addPair`
/// makes them, each with the other's CIGAR in its `MC` tag.
fn mates(start1: usize, cigar1: &str, start2: usize, cigar2: &str) -> Vec<RecordBuf> {
    SamBuilder::new().add_pair(Pair {
        name: Some("q".into()),
        cigar1: Some(cigar1.into()),
        cigar2: Some(cigar2.into()),
        ..Pair::at(start1, start2)
    })
}

/// The template-end distance of the forward or the reverse read at each 0-based position.
fn ends(records: &[RecordBuf], reverse: bool, positions: &[usize]) -> Vec<Option<usize>> {
    let mut builder = piled(records.to_vec());
    positions
        .iter()
        .map(|&position| {
            let pileup = builder.pileup("chr1", position).unwrap();
            let entry = pileup
                .iter()
                .find(|entry| entry.is_reverse() == reverse)
                .unwrap();
            entry.template_end_distance().unwrap()
        })
        .collect()
}

/// fgbio's `SamRecordTest`: "SamRecord.mateUnclippedStart/End should return None if no mate cigar
/// set". A missing `MC` tag is an error naming the read here, for either read of the pair.
#[test]
fn test_a_template_end_without_a_mate_cigar_is_an_error_for_either_read() {
    let records = pair("q", 10, 50, Strand::Plus, Strand::Minus)
        .into_iter()
        .map(SamBuilder::without_mate_cigar);
    let mut builder = piled(records);
    let pileup = builder.pileup("chr1", 50).unwrap();
    assert_eq!(pileup.len(), 2);
    for entry in pileup.iter() {
        let refused = entry.template_end_distance();
        assert!(
            matches!(&refused, Err(Error::MissingMateCigar { name }) if name == "q"),
            "{refused:?}"
        );
    }
}

/// fgbio's `SamRecordTest` fixtures for the mate's unclipped start and end, counted in template
/// bases: a mate's soft-clipped bases count and its hard-clipped bases, absent from it, do not.
#[test]
fn test_a_template_end_counts_soft_clips_and_not_hard_clips() {
    let reads = mates(10, "5S45M10H", 50, "10S40M5H");
    assert_eq!(ends(&reads, false, &[50]), [Some(88 - 50)]);
    assert_eq!(ends(&reads, true, &[50]), [Some(5 + 50 - 9)]);
    let reads = mates(10, "10H5S45M1S10H", 50, "10H10S40M1S5H");
    assert_eq!(ends(&reads, false, &[50]), [Some(88 + 1 - 50)]);
    assert_eq!(ends(&reads, true, &[50]), [Some(5 + 50 - 9)]);
}

/// fgbio #1172: "clip back to the mate's soft-clipped end when neither read contains an indel".
/// The forward read's last 40 bases and the reverse read's first 40 lie past the template's ends.
#[test]
fn test_a_template_end_without_indels_reaches_the_mates_soft_clipped_end() {
    let reads = mates(100, "100M", 60, "90M10S");
    assert_eq!(
        ends(&reads, false, &[157, 158, 159]),
        [Some(1), Some(0), None]
    );
    assert_eq!(ends(&reads, true, &[98, 99]), [None, Some(0)]);
}

/// fgbio #1172: "clip using query distances when a deletion falls at the mate's un-soft-clipped
/// end". Past the last shared position, 222, the forward read has a base, a deleted position, and
/// three bases, and the mate two bases, so the first of the three is the template's last base.
#[test]
fn test_a_template_end_counts_bases_past_a_deletion_at_the_mates_soft_clipped_end() {
    let reads = mates(101, "2S124M1D3M", 100, "3S124M2S");
    assert_eq!(
        ends(&reads, false, &[223, 224, 225, 226]),
        [Some(1), Some(1), Some(0), None]
    );
    assert_eq!(ends(&reads, true, &[99, 100]), [Some(1), Some(2)]);
}

/// fgbio #1172: "clip the mate even when the record's deletion falls at the mate's un-soft-clipped
/// end". The reverse read's first two aligned bases lie before the forward read's 5′ end.
#[test]
fn test_a_template_end_is_found_for_both_reads_when_a_deletion_meets_the_mates_end() {
    let reads = mates(101, "2S124M1D3M", 97, "115M14S");
    assert_eq!(
        ends(&reads, false, &[223, 225, 226]),
        [Some(1), Some(0), None]
    );
    assert_eq!(
        ends(&reads, true, &[96, 97, 98, 99]),
        [None, None, Some(0), Some(1)]
    );
}

/// fgbio #1172: "clip using query distances when the record contains an insertion before the
/// mate's end". The 10 inserted bases bring the template's end 10 reference positions closer.
#[test]
fn test_a_template_end_counts_inserted_bases() {
    let reads = mates(100, "70M10I23M47S", 100, "50S70M30S");
    assert_eq!(
        ends(&reads, false, &[168, 169, 188, 189]),
        [Some(30), Some(19), Some(0), None]
    );
}

/// fgbio #1172: "clip the read-through of a pair whose alignments share no reference position".
/// Moving the mate one base right lengthens the template by one, whether or not the reads share a
/// reference position.
#[test]
fn test_a_template_end_is_continuous_as_the_reads_stop_sharing_a_position() {
    let reads = mates(1000, "20M80S", 1019, "80S20M");
    assert_eq!(ends(&reads, false, &[1018]), [Some(19)]);
    assert_eq!(ends(&reads, true, &[1018]), [Some(19)]);
    let reads = mates(1000, "20M80S", 1020, "80S20M");
    assert_eq!(ends(&reads, false, &[1018]), [Some(20)]);
    assert_eq!(ends(&reads, true, &[1019]), [Some(20)]);
}

/// fgbio #1172: "not clip a record lying entirely past the mate alignment it is given". A forward
/// read whose mate lies wholly before it has no template end at any of its positions.
#[test]
fn test_a_template_end_is_absent_for_a_read_past_its_mates_alignment() {
    let records = mates(301, "15M", 101, "10M5S");
    assert_eq!(ends(&records, false, &[300, 307, 314]), [None, None, None]);
}

#[test]
fn test_a_mate_cigar_must_span_a_base_with_operators_bam_allows() {
    let mate = |text: &str| {
        parse_cigar(text.as_bytes())
            .map(|ops| Alignment::new(0, ops))
            .filter(Alignment::spans_reference)
    };
    for valid in ["10M2I3D4N1=1X5S2H1P", "+3M", "268435455M268435455N"] {
        assert!(mate(valid).is_some(), "{valid}");
    }
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
        assert_eq!(mate(invalid), None, "{invalid}");
    }
}

/// The template ends of any alignment record, here a `RecordBuf` holding fgbio #1172's deletion
/// case, as the public functions give them.
#[test]
fn test_the_template_ends_of_a_record_buf() -> crate::Result<()> {
    let mut builder = SamBuilder::new();
    let record = builder.add_pair(Pair {
        cigar1: Some("2S124M1D3M".into()),
        cigar2: Some("3S124M2S".into()),
        ..Pair::at(101, 100)
    })[0]
        .clone();
    let header = builder.header();
    let ends: Vec<_> = [223, 224, 225, 226]
        .into_iter()
        .map(|position| crate::template_end_distance(&record, header, position).unwrap())
        .collect();
    assert_eq!(ends, [Some(1), Some(1), Some(0), None]);
    let five_prime =
        |offset| -> crate::Result<Option<usize>> { crate::five_prime_distance(&record, offset) };
    assert_eq!(
        (five_prime(0)?, five_prime(128)?, five_prime(129)?),
        (Some(0), Some(128), None)
    );
    Ok(())
}
