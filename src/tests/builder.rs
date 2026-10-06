use std::collections::HashMap;
use std::io;
use std::sync::{Arc, Mutex};

use noodles::bam;
use noodles::sam::alignment::record::Flags;

use super::{
    Entry, HEADER, Read, bam_bytes, bases, builder, entries, entry, name, names, raw_record, read,
    records, source, unmapped,
};
use crate::{EntryKind, Error, PileupEntry, RecordSource, StreamingPileupBuilder};

type Log = Arc<Mutex<Vec<String>>>;

fn log() -> Log {
    Arc::new(Mutex::new(Vec::new()))
}

fn logged(log: &Log) -> Vec<String> {
    log.lock().unwrap().clone()
}

/// A tap that notes the name of each record it is given.
fn tap_into(log: &Log) -> impl FnMut(bam::Record) -> io::Result<()> + Send + 'static {
    let log = Arc::clone(log);
    move |record| {
        log.lock().unwrap().push(name(&record));
        Ok(())
    }
}

/// A source that notes the name of each record as it is read.
struct Counted {
    records: std::vec::IntoIter<bam::Record>,
    pulled: Log,
}

impl RecordSource for Counted {
    fn read_record(&mut self, record: &mut bam::Record) -> io::Result<bool> {
        match self.records.next() {
            Some(next) => {
                self.pulled.lock().unwrap().push(name(&next));
                *record = next;
                Ok(true)
            }
            None => Ok(false),
        }
    }
}

fn read_names(reads: &[Read]) -> Vec<String> {
    let (_, records) = records(HEADER, reads);
    records.iter().map(name).collect()
}

#[test]
fn test_builder_refuses_records_not_declared_coordinate_sorted() {
    let header = "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:100\n";
    let (header, queryname) = records(header, &[read("r", 10, "4M", "ACGT")]);
    let refused = StreamingPileupBuilder::new(source(queryname), &header).err();
    assert!(
        matches!(refused, Some(Error::NotCoordinateSorted { found: Some(order) }) if order == "queryname")
    );
    let header = "@SQ\tSN:chr1\tLN:100\n";
    let (header, unsorted) = records(header, &[]);
    let refused = StreamingPileupBuilder::new(source(unsorted), &header).err();
    assert!(matches!(
        refused,
        Some(Error::NotCoordinateSorted { found: None })
    ));
}

#[test]
fn test_builder_refuses_records_out_of_order() {
    let mut builder = builder(&[read("b", 20, "4M", "ACGT"), read("a", 10, "4M", "ACGT")]);
    let error = builder.pileup("chr1", 30).err().unwrap();
    assert_eq!(
        error.to_string(),
        "records are out of coordinate order at a"
    );
}

#[test]
fn test_builder_is_forward_only() {
    let mut builder = builder(&[read("r", 100, "4M", "ACGT")]);
    builder.pileup("chr2", 50).unwrap();
    let error = builder.pileup("chr1", 100).err().unwrap();
    assert_eq!(
        error.to_string(),
        "attempted to advance to chr1:100 from chr2:50"
    );
    let error = builder.pileup("chr2", 20).err().unwrap();
    assert_eq!(
        error.to_string(),
        "attempted to advance to chr2:20 from chr2:50"
    );
}

#[test]
fn test_builder_refuses_unknown_contigs_and_spans_that_end_before_they_start() {
    let mut builder = builder(&[read("r", 100, "4M", "ACGT")]);
    assert_eq!(
        builder.pileup("chr3", 1).err().unwrap().to_string(),
        "contig chr3 is not in the header"
    );
    assert!(matches!(
        builder.pileup_at(2, 1),
        Err(Error::UnknownReferenceSequenceId(2))
    ));
    assert!(matches!(
        builder.columns("chr1", 2, 1).err(),
        Some(Error::InvalidSpan { start: 2, end: 1 })
    ));
}

#[test]
fn test_builder_refuses_pileups_once_closed() {
    let mut builder = builder(&[read("r", 10, "4M", "ACGT")]);
    builder.pileup("chr1", 11).unwrap();
    builder.close().unwrap();
    for position in [11, 12] {
        assert!(matches!(
            builder.pileup("chr1", position),
            Err(Error::Closed)
        ));
    }
    let mut columns = builder.columns("chr1", 12, 14).unwrap();
    assert!(matches!(columns.next_pileup(), Some(Err(Error::Closed))));
    builder.close().unwrap();
}

#[test]
fn test_builder_returns_the_same_pileup_for_a_repeated_position() {
    let pulled = log();
    let (header, records) = records(
        HEADER,
        &[read("r", 100, "4M", "ACGT"), read("s", 102, "4M", "ACGT")],
    );
    let counted = Counted {
        records: records.into_iter(),
        pulled: Arc::clone(&pulled),
    };
    let mut builder = StreamingPileupBuilder::new(counted, &header).unwrap();
    let first = entries(&builder.pileup("chr1", 101).unwrap());
    let read_so_far = logged(&pulled);
    let again = entries(&builder.pileup("chr1", 101).unwrap());
    assert_eq!(first, [entry("r", "base", Some(1), Some(1), None)]);
    assert_eq!((again, logged(&pulled)), (first, read_so_far));
}

#[test]
fn test_builder_checks_the_header_of_an_empty_input() {
    let mut builder = builder(&[]);
    assert!(builder.pileup("chr2", 5).unwrap().is_empty());
    assert!(matches!(
        builder.pileup("chr3", 1),
        Err(Error::UnknownContig(_))
    ));
    let error = builder.pileup("chr1", 5).err().unwrap();
    assert_eq!(
        error.to_string(),
        "attempted to advance to chr1:5 from chr2:5"
    );
}

#[test]
fn test_builder_piles_up_a_deletion() {
    let mut builder = builder(&[read("r", 10, "2M2D2M", "ACGT").quals(&[30, 31, 32, 33])]);
    assert_eq!(
        entries(&builder.pileup("chr1", 11).unwrap()),
        [entry("r", "base", Some(1), Some(1), None)]
    );
    let pileup = builder.pileup("chr1", 12).unwrap();
    assert_eq!(
        entries(&pileup),
        [entry("r", "deletion", None, Some(2), None)]
    );
    let deletion = pileup.get(0).unwrap();
    assert!(deletion.is_deletion() && deletion.quality() == Some(32) && deletion.base().is_none());
    assert_eq!(
        entries(&builder.pileup("chr1", 13).unwrap()),
        [entry("r", "deletion", None, Some(2), None)]
    );
    assert_eq!(
        entries(&builder.pileup("chr1", 14).unwrap()),
        [entry("r", "base", Some(2), Some(2), None)]
    );
}

#[test]
fn test_builder_piles_up_a_deletion_no_base_follows() {
    let reads = [read("r", 10, "3M2D", "ACG"), read("s", 10, "4M", "ACGT")];
    let mut builder = builder(&reads).min_base_quality(0);
    let pileup = builder.pileup("chr1", 13).unwrap();
    assert_eq!(
        entries(&pileup),
        [
            entry("r", "deletion", None, None, None),
            entry("s", "base", Some(3), Some(3), None)
        ]
    );
    assert_eq!(pileup.get(0).unwrap().quality(), None);
    assert_eq!((pileup.unfiltered_depth(), pileup.filtered_depth()), (2, 1));
    assert_eq!(
        entries(&builder.pileup("chr1", 14).unwrap()),
        [entry("r", "deletion", None, None, None)]
    );
    assert_eq!(
        entries(&builder.pileup("chr1", 15).unwrap()),
        Vec::<Entry>::new()
    );
}

#[test]
fn test_builder_piles_up_a_read_that_opens_with_a_deletion() {
    let cases: [(&str, &str, Vec<Entry>); 3] = [
        (
            "1D3M",
            "ACG",
            vec![entry("r", "deletion", None, Some(0), None)],
        ),
        (
            "3D4M",
            "ACGT",
            vec![entry("r", "deletion", None, Some(0), None)],
        ),
        (
            "1D1I3M",
            "TACG",
            vec![
                entry("r", "deletion", None, Some(0), None),
                entry("r", "insertion", None, None, Some("T")),
            ],
        ),
    ];
    for (cigar, read_bases, expected) in cases {
        let mut quals = vec![40; read_bases.len()];
        quals[0] = 25;
        let mut builder = builder(&[read("r", 10, cigar, read_bases).quals(&quals)]);
        assert!(builder.pileup("chr1", 9).unwrap().is_empty());
        let pileup = builder.pileup("chr1", 10).unwrap();
        assert_eq!(entries(&pileup), expected, "{cigar}");
        assert_eq!(pileup.get(0).unwrap().quality(), Some(25));
    }
}

#[test]
fn test_builder_piles_up_an_opening_insertion_after_a_hard_clip() {
    let mut builder = builder(&[read("r", 10, "4H1I3M", "TACG")]);
    assert_eq!(
        entries(&builder.pileup("chr1", 9).unwrap()),
        [entry("r", "insertion", None, None, Some("T"))]
    );
    assert_eq!(
        entries(&builder.pileup("chr1", 10).unwrap()),
        [entry("r", "base", Some(1), Some(1), None)]
    );
}

#[test]
fn test_builder_counts_a_column_of_every_kind_of_entry() {
    let reads = [
        read("base", 10, "4M", "ACGT"),
        read("lowdel", 10, "2M1D2M", "ACTA").quals(&[40, 40, 5, 40]),
        read("del", 10, "2M1D2M", "ACTA"),
        read("closes", 10, "3M1I1M", "ACGTA"),
        read("opens", 13, "1I3M", "TACG"),
    ];
    let mut builder = builder(&reads);
    let pileup = builder.pileup("chr1", 12).unwrap();
    assert_eq!(
        entries(&pileup),
        [
            entry("base", "base", Some(2), Some(2), None),
            entry("lowdel", "deletion", None, Some(2), None),
            entry("del", "deletion", None, Some(2), None),
            entry("closes", "base", Some(2), Some(2), None),
            entry("closes", "insertion", None, None, Some("T")),
            entry("opens", "insertion", None, None, Some("T")),
        ]
    );
    assert_eq!(
        (
            pileup.len(),
            pileup.unfiltered_depth(),
            pileup.filtered_depth()
        ),
        (6, 4, 3)
    );
    assert_eq!(bases(&pileup), "GG");
}

#[test]
fn test_builder_floor_leaves_bases_under_it_out_of_the_views_only() {
    let reads: Vec<Read> = [19, 20, 21]
        .map(|quality| {
            read(&format!("q{quality}"), 100, "50M", &"A".repeat(50)).quals(&[quality; 50])
        })
        .into();
    let mut builder = builder(&reads).min_base_quality(20);
    let pileup = builder.pileup("chr1", 104).unwrap();
    assert_eq!(
        (
            pileup.filtered_depth(),
            pileup.qualities().collect::<Vec<_>>()
        ),
        (2, vec![20, 21])
    );
    assert_eq!(bases(&pileup), "AA");
    assert_eq!(names(&pileup), ["q19", "q20", "q21"]);
    assert_eq!(pileup.unfiltered_depth(), 3);
}

#[test]
fn test_builder_piles_up_a_crowd_of_reads_at_one_start() {
    let mut reads = Vec::new();
    for (base, count) in [('A', 5), ('C', 4), ('G', 3), ('T', 2), ('N', 1)] {
        for index in 0..count {
            reads.push(read(
                &format!("{base}{index}"),
                5,
                "10M",
                &base.to_string().repeat(10),
            ));
        }
    }
    let mut builder = builder(&reads);
    let pileup = builder.pileup("chr1", 5).unwrap();
    let mut counts: HashMap<char, usize> = HashMap::new();
    for base in bases(&pileup).chars() {
        *counts.entry(base).or_default() += 1;
    }
    assert_eq!(
        counts,
        HashMap::from([('A', 5), ('C', 4), ('G', 3), ('T', 2), ('N', 1)])
    );
    assert_eq!(pileup.iter().filter(PileupEntry::is_no_call).count(), 1);
    assert_eq!(pileup.unfiltered_depth(), 15);
    assert_eq!(names(&pileup), read_names(&reads));
}

#[test]
fn test_builder_piles_up_read_through_pairs_and_reverse_reads_as_aligned() {
    let reads = [
        read("pair", 99, "10M", "ACGTACGTAC")
            .flag(83)
            .quals(&[35; 10]),
        read("pair", 100, "10M", "CGTACGTACG")
            .flag(163)
            .quals(&[35; 10]),
    ];
    let mut builder = builder(&reads);
    let mut columns = builder.columns("chr1", 98, 111).unwrap();
    let mut depths = Vec::new();
    while let Some(pileup) = columns.next_pileup() {
        depths.push(pileup.unwrap().unfiltered_depth());
    }
    assert_eq!(depths, [0, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 0]);
    let mut builder = super::builder(&reads);
    let pileup = builder.pileup("chr1", 103).unwrap();
    let reverse = pileup.get(0).unwrap();
    assert!(reverse.is_reverse());
    assert_eq!(
        (reverse.query_position(), reverse.base(), reverse.quality()),
        (Some(4), Some(b'A'), Some(35))
    );
}

#[test]
fn test_builder_piles_up_insertions_at_read_starts_and_ends() {
    let reads = [
        read("opens", 10, "2I4M", "TTACGT").quals(&[10, 11, 30, 30, 30, 30]),
        read("closes", 10, "4M2I", "ACGTCC"),
        read("clipped", 10, "1S1I4M", "GTACGT"),
    ];
    let mut builder = builder(&reads);
    let at_nine = builder.pileup("chr1", 9).unwrap();
    assert_eq!(
        entries(&at_nine),
        [
            entry("opens", "insertion", None, None, Some("TT")),
            entry("clipped", "insertion", None, None, Some("T"))
        ]
    );
    let opens = at_nine.get(0).unwrap();
    assert!(opens.is_insertion());
    assert_eq!(
        opens.inserted_qualities().map(Iterator::collect::<Vec<_>>),
        Some(vec![10, 11])
    );
    let offsets: Vec<Option<usize>> = at_nine
        .iter()
        .map(|entry| entry.insertion_offset())
        .collect();
    assert_eq!(offsets, [Some(0), Some(1)]);
    assert_eq!(at_nine.unfiltered_depth(), 0);
    assert_eq!(
        entries(&builder.pileup("chr1", 10).unwrap()),
        [
            entry("opens", "base", Some(2), Some(2), None),
            entry("closes", "base", Some(0), Some(0), None),
            entry("clipped", "base", Some(2), Some(2), None),
        ]
    );
    assert_eq!(
        entries(&builder.pileup("chr1", 13).unwrap()),
        [
            entry("opens", "base", Some(5), Some(5), None),
            entry("closes", "base", Some(3), Some(3), None),
            entry("closes", "insertion", None, None, Some("CC")),
            entry("clipped", "base", Some(5), Some(5), None),
        ]
    );
    assert!(builder.pileup("chr1", 14).unwrap().is_empty());
}

#[test]
fn test_reads_with_no_stored_qualities_pass_every_floor() {
    let mut builder =
        builder(&[read("r", 10, "2M1D1M2I", "ACGTT").no_quals()]).min_base_quality(60);
    let at_base = builder.pileup("chr1", 11).unwrap();
    assert_eq!(bases(&at_base), "C");
    assert_eq!(
        (
            at_base.qualities().collect::<Vec<_>>(),
            at_base.filtered_depth()
        ),
        (vec![255], 1)
    );
    let at_deletion = builder.pileup("chr1", 12).unwrap();
    assert_eq!(
        (
            at_deletion.get(0).unwrap().quality(),
            at_deletion.filtered_depth()
        ),
        (Some(255), 1)
    );
    let at_insertion = builder.pileup("chr1", 13).unwrap();
    assert_eq!(
        entries(&at_insertion),
        [
            entry("r", "base", Some(2), Some(2), None),
            entry("r", "insertion", None, None, Some("TT"))
        ]
    );
    let inserted = at_insertion
        .get(1)
        .unwrap()
        .inserted_qualities()
        .map(Iterator::collect::<Vec<_>>);
    assert_eq!(inserted, Some(vec![255, 255]));
}

#[test]
fn test_reads_with_no_stored_bases_hold_no_base_or_quality() {
    let mut builder = builder(&[read("r", 10, "2M1I1M", "*")]);
    let pileup = builder.pileup("chr1", 11).unwrap();
    assert_eq!(
        entries(&pileup),
        [
            entry("r", "base", Some(1), Some(1), None),
            entry("r", "insertion", None, None, None)
        ]
    );
    let held: Vec<(Option<u8>, Option<u8>)> = pileup
        .iter()
        .map(|entry| (entry.base(), entry.quality()))
        .collect();
    assert_eq!(held, [(None, None), (None, None)]);
    assert!(pileup.get(1).unwrap().inserted_qualities().is_none());
    assert_eq!((pileup.unfiltered_depth(), pileup.filtered_depth()), (1, 0));
    assert_eq!(
        (bases(&pileup), pileup.qualities().count()),
        (String::new(), 0)
    );
}

#[test]
fn test_builder_leaves_out_reads_with_no_reference_consuming_operator() {
    let reads = [
        read("clipped", 10, "2S2I", "ACGT"),
        read("inserted", 10, "4I", "ACGT"),
    ];
    let evicted = log();
    let mut builder = builder(&reads).tap(tap_into(&evicted));
    let mut columns = builder.columns("chr1", 8, 12).unwrap();
    while let Some(pileup) = columns.next_pileup() {
        assert!(pileup.unwrap().is_empty());
    }
    builder.close().unwrap();
    assert_eq!(logged(&evicted), ["clipped", "inserted"]);
}

#[test]
fn test_builder_skips_soft_and_hard_clips() {
    let mut builder = builder(&[read("r", 10, "5H2S3M1S", "TTACGA")]);
    assert!(builder.pileup("chr1", 9).unwrap().is_empty());
    let pileup = builder.pileup("chr1", 10).unwrap();
    assert_eq!(
        entries(&pileup),
        [entry("r", "base", Some(2), Some(2), None)]
    );
    assert_eq!(bases(&pileup), "A");
    assert_eq!(
        entries(&builder.pileup("chr1", 12).unwrap()),
        [entry("r", "base", Some(4), Some(4), None)]
    );
    assert!(builder.pileup("chr1", 13).unwrap().is_empty());
}

#[test]
fn test_builder_piles_up_reference_skips() {
    let mut builder = builder(&[read("r", 10, "2M3N2M", "ACGT")]).min_base_quality(0);
    let mut columns = builder.columns("chr1", 10, 18).unwrap();
    let mut swept = Vec::new();
    while let Some(pileup) = columns.next_pileup() {
        let pileup = pileup.unwrap();
        let skip = pileup.get(0).filter(PileupEntry::is_skip);
        if let Some(skip) = skip {
            assert!(!skip.is_deletion() && !skip.is_insertion());
            assert_eq!((skip.base(), skip.quality()), (None, None));
            assert_eq!(
                (bases(&pileup), pileup.qualities().count()),
                (String::new(), 0)
            );
        }
        swept.push((
            entries(&pileup),
            pileup.unfiltered_depth(),
            pileup.filtered_depth(),
        ));
    }
    let skip = || vec![entry("r", "skip", None, None, None)];
    assert_eq!(
        swept,
        [
            (vec![entry("r", "base", Some(0), Some(0), None)], 1, 1),
            (vec![entry("r", "base", Some(1), Some(1), None)], 1, 1),
            (skip(), 1, 0),
            (skip(), 1, 0),
            (skip(), 1, 0),
            (vec![entry("r", "base", Some(2), Some(2), None)], 1, 1),
            (vec![entry("r", "base", Some(3), Some(3), None)], 1, 1),
            (vec![], 0, 0),
        ]
    );
}

#[test]
fn test_builder_piles_up_an_insertion_after_a_reference_skip() {
    let mut builder = builder(&[read("r", 10, "2M2N1I2M", "ACTGT")]);
    assert_eq!(
        entries(&builder.pileup("chr1", 13).unwrap()),
        [
            entry("r", "skip", None, None, None),
            entry("r", "insertion", None, None, Some("T"))
        ]
    );
    assert_eq!(
        entries(&builder.pileup("chr1", 14).unwrap()),
        [entry("r", "base", Some(3), Some(3), None)]
    );
}

#[test]
fn test_builder_piles_up_both_overlapping_mates() {
    let reads = [
        read("pair", 10, "6M", "ACGTAC").flag(99),
        read("pair", 13, "6M", "TACGTA").flag(147),
    ];
    let mut builder = builder(&reads);
    let pileup = builder.pileup("chr1", 14).unwrap();
    let held: Vec<(u16, Option<usize>)> = pileup
        .iter()
        .map(|entry| (entry.flags().bits(), entry.query_position()))
        .collect();
    assert_eq!(held, [(99, Some(4)), (147, Some(1))]);
    assert_eq!(bases(&pileup), "AA");
}

#[test]
fn test_builder_moves_across_contigs() {
    let reads = [
        read("one", 10, "4M", "ACGT"),
        read("two", 10, "4M", "TTTT").contig("chr2"),
    ];
    let evicted = log();
    let mut builder = builder(&reads).tap(tap_into(&evicted));
    assert_eq!(bases(&builder.pileup("chr1", 11).unwrap()), "C");
    assert_eq!(bases(&builder.pileup("chr2", 11).unwrap()), "T");
    assert_eq!(logged(&evicted), ["one"]);
    builder.close().unwrap();
    assert_eq!(logged(&evicted), ["one", "two"]);
}

#[test]
fn test_builder_floors_bases_at_13_and_leaves_out_qc_fail_reads_by_default() {
    let reads = [
        read("q12", 10, "4M", "ACGT").quals(&[12; 4]),
        read("q13", 10, "4M", "ACGT").quals(&[13; 4]),
        read("qcfail", 10, "4M", "ACGT").flag(512),
    ];
    let mut builder = builder(&reads);
    let pileup = builder.pileup("chr1", 10).unwrap();
    assert_eq!(names(&pileup), ["q12", "q13"]);
    assert_eq!(
        (
            pileup.filtered_depth(),
            pileup.qualities().collect::<Vec<_>>()
        ),
        (1, vec![13])
    );
}

#[test]
fn test_builder_filters_reads() {
    type Options = fn(
        StreamingPileupBuilder<'static, super::Source>,
    ) -> StreamingPileupBuilder<'static, super::Source>;
    let cases: [(u16, u8, Options, bool); 14] = [
        (0, 5, |b| b.min_mapping_quality(20), false),
        (0, 20, |b| b.min_mapping_quality(20), true),
        (256, 60, |b| b, false),
        (
            256,
            60,
            |b| b.exclude_flags(Flags::from_bits_retain(0xE00)),
            true,
        ),
        (2048, 60, |b| b, false),
        (
            2048,
            60,
            |b| b.exclude_flags(Flags::from_bits_retain(0x700)),
            true,
        ),
        (1024, 60, |b| b, false),
        (
            1024,
            60,
            |b| b.exclude_flags(Flags::from_bits_retain(0xB00)),
            true,
        ),
        (512, 60, |b| b, false),
        (
            512,
            60,
            |b| b.exclude_flags(Flags::from_bits_retain(0xD00)),
            true,
        ),
        (
            512,
            60,
            |b| b.exclude_flags(Flags::from_bits_retain(0xE00)),
            false,
        ),
        (
            16,
            60,
            |b| b.exclude_flags(Flags::from_bits_retain(0x10)),
            false,
        ),
        (1, 60, |b| b.proper_pairs_only(true), false),
        (3, 60, |b| b.proper_pairs_only(true), true),
    ];
    for (flag, mapq, options, kept) in cases {
        let evicted = log();
        let mut builder = options(builder(&[read("r", 10, "4M", "ACGT")
            .flag(flag)
            .mapq(mapq)]))
        .tap(tap_into(&evicted));
        assert_eq!(
            builder.pileup("chr1", 10).unwrap().unfiltered_depth(),
            usize::from(kept),
            "{flag} {mapq}"
        );
        builder.close().unwrap();
        assert_eq!(logged(&evicted), ["r"]);
    }
}

#[test]
fn test_a_read_with_no_mapping_quality_passes_every_floor() {
    let mut builder = builder(&[read("r", 10, "4M", "ACGT").mapq(255)]).min_mapping_quality(60);
    assert_eq!(builder.pileup("chr1", 10).unwrap().unfiltered_depth(), 1);
}

#[test]
fn test_builder_asks_a_read_filter_after_its_own_filters() {
    let reads = [
        read("q1", 100, "50M", &"A".repeat(50)),
        read("x2", 104, "50M", &"A".repeat(50)),
        read("q3", 108, "50M", &"A".repeat(50)).flag(1024),
        read("q4", 112, "50M", &"A".repeat(50)),
        unmapped("q5"),
    ];
    let asked = log();
    let evicted = log();
    let asking = Arc::clone(&asked);
    let mut builder = builder(&reads)
        .read_filter(move |record| {
            let read_name = name(record);
            asking.lock().unwrap().push(read_name.clone());
            !read_name.starts_with('x')
        })
        .tap(tap_into(&evicted));
    assert_eq!(names(&builder.pileup("chr1", 115).unwrap()), ["q1", "q4"]);
    builder.close().unwrap();
    assert_eq!(logged(&asked), ["q1", "x2", "q4"]);
    assert_eq!(logged(&evicted), read_names(&reads));
}

#[test]
fn test_builder_taps_every_read_once_in_input_order() {
    let reads = [
        read("long", 100, "50M", &"A".repeat(50)),
        read("short", 100, "40M", &"A".repeat(40)),
        read("filtered", 110, "10M", &"A".repeat(10)).flag(1024),
        read("later", 200, "50M", &"A".repeat(50)),
        read("other", 10, "4M", "ACGT").contig("chr2"),
        unmapped("unplaced"),
    ];
    let evicted = log();
    let mut builder = builder(&reads).tap(tap_into(&evicted));
    builder.pileup("chr1", 100).unwrap();
    builder.pileup("chr1", 145).unwrap();
    assert_eq!(logged(&evicted), Vec::<String>::new());
    builder.pileup("chr1", 200).unwrap();
    assert_eq!(logged(&evicted), ["long", "short", "filtered"]);
    builder.close().unwrap();
    assert_eq!(logged(&evicted), read_names(&reads));
}

#[test]
fn test_builder_taps_a_read_once_it_has_passed_the_read_and_every_one_before() {
    let reads = [
        read("first", 10, "4M", "ACGT"),
        read("unmapped", 10, "*", "ACGT").flag(4),
        read("second", 12, "4M", "ACGT"),
    ];
    let evicted = log();
    let mut builder = builder(&reads).tap(tap_into(&evicted));
    assert_eq!(
        names(&builder.pileup("chr1", 13).unwrap()),
        ["first", "second"]
    );
    assert_eq!(logged(&evicted), Vec::<String>::new());
    builder.pileup("chr1", 14).unwrap();
    assert_eq!(logged(&evicted), ["first", "unmapped"]);
}

#[test]
fn test_builder_taps_every_read_when_dropped_after_an_error() {
    let reads: Vec<Read> = [10, 20, 30]
        .map(|start| read(&format!("r{start}"), start, "4M", "ACGT"))
        .into();
    let evicted = log();
    let failing = || -> Result<(), Error> {
        let mut builder = builder(&reads).tap(tap_into(&evicted));
        builder.pileup("chr1", 10)?;
        builder.pileup("chr3", 10)?;
        Ok(())
    };
    assert!(failing().is_err());
    assert_eq!(logged(&evicted), ["r10", "r20", "r30"]);
}

#[test]
fn test_closing_reads_the_rest_of_the_input_only_for_a_tap() {
    for tapped in [false, true] {
        let reads: Vec<Read> = (1..6)
            .map(|index| read(&format!("r{}", index * 10), index * 10, "4M", "ACGT"))
            .collect();
        let (header, records) = records(HEADER, &reads);
        let pulled = log();
        let evicted = log();
        let counted = Counted {
            records: records.into_iter(),
            pulled: Arc::clone(&pulled),
        };
        let mut builder = StreamingPileupBuilder::new(counted, &header).unwrap();
        if tapped {
            builder = builder.tap(tap_into(&evicted));
        }
        builder.pileup("chr1", 10).unwrap();
        assert_eq!(logged(&pulled), ["r10", "r20"]);
        drop(builder);
        let all = read_names(&reads);
        assert_eq!(
            logged(&pulled),
            if tapped {
                all.clone()
            } else {
                vec!["r10".into(), "r20".into()]
            }
        );
        assert_eq!(logged(&evicted), if tapped { all } else { vec![] });
    }
}

#[test]
fn test_builder_holds_passed_reads_only_for_a_tap() {
    let reads = [
        read("long", 100, "100M", &"A".repeat(100)),
        read("short", 100, "10M", &"A".repeat(10)),
        read("filtered", 105, "10M", &"A".repeat(10)).flag(1024),
        read("later", 300, "4M", "ACGT"),
    ];
    for tapped in [false, true] {
        let mut builder = builder(&reads);
        if tapped {
            builder = builder.tap(|_| Ok(()));
        }
        builder.pileup("chr1", 120).unwrap();
        let mut held: Vec<String> = builder
            .active
            .iter()
            .chain(builder.waiting.iter())
            .chain(builder.next.iter())
            .map(|&index| name(&builder.slots[index as usize].record))
            .collect();
        held.sort();
        held.dedup();
        let expected = if tapped {
            vec!["filtered", "later", "long", "short"]
        } else {
            vec!["later", "long"]
        };
        assert_eq!(held, expected);
    }
}

#[test]
fn test_builder_hands_records_to_a_tap_that_writes_them() {
    let reads = [
        read("one", 100, "4M", "ACGT"),
        read("two", 200, "4M", "ACGT"),
    ];
    let (header, records) = records(HEADER, &reads);
    let written = Arc::new(Mutex::new(Vec::new()));
    let sink = Arc::clone(&written);
    let mut builder = StreamingPileupBuilder::new(source(records), &header)
        .unwrap()
        .tap(move |record| {
            sink.lock().unwrap().push(record);
            Ok(())
        });
    for position in [100, 200] {
        assert_eq!(builder.pileup("chr1", position).unwrap().len(), 1);
    }
    builder.close().unwrap();
    let mut writer = bam::io::Writer::new(Vec::new());
    writer.write_header(&header).unwrap();
    for record in written.lock().unwrap().iter() {
        writer.write_record(&header, record).unwrap();
    }
    let bytes = writer.into_inner().finish().unwrap();
    let mut reader = bam::io::Reader::new(&bytes[..]);
    reader.read_header().unwrap();
    let round_trip: Vec<String> = reader
        .records()
        .map(|record| name(&record.unwrap()))
        .collect();
    assert_eq!(round_trip, ["one", "two"]);
}

#[test]
fn test_a_tap_error_is_returned() {
    let reads = [read("one", 10, "4M", "ACGT"), read("two", 20, "4M", "ACGT")];
    let mut builder = builder(&reads).tap(|_| Err(io::Error::other("full")));
    assert!(builder.pileup("chr1", 10).is_ok());
    assert!(matches!(builder.pileup("chr1", 20), Err(Error::Io(_))));
}

#[test]
fn test_columns_agree_with_pileups_at_every_position() {
    let reads = [
        read("a", 5, "3S10M2D5M", "GGGACGTACGTACTTTTT"),
        read("b", 7, "4M2I6M", "ACGTGGACGTAC"),
        read("c", 9, "2I8M", "TTACGTACGT"),
        read("d", 10, "3M4N3M", "ACGTAC"),
        read("e", 12, "8M", "ACGTACGT").flag(1024),
        read("f", 14, "3M2D", "ACG"),
        read("g", 30, "5M", "ACGTA").contig("chr2"),
    ];
    let spans = [("chr1", 0, 30), ("chr2", 25, 40)];
    let mut swept = Vec::new();
    for (contig, start, end) in spans {
        let mut builder = builder(&reads);
        let mut columns = builder.columns(contig, start, end).unwrap();
        while let Some(pileup) = columns.next_pileup() {
            swept.push(entries(&pileup.unwrap()));
        }
    }
    let mut one_by_one = Vec::new();
    for (contig, start, end) in spans {
        for position in start..end {
            one_by_one.push(entries(&builder(&reads).pileup(contig, position).unwrap()));
        }
    }
    let kept: Vec<Read> = reads
        .iter()
        .filter(|read| read.flag & 1024 == 0)
        .cloned()
        .collect();
    let mut unfiltered = Vec::new();
    for (contig, start, end) in spans {
        for position in start..end {
            unfiltered.push(entries(
                &super::unfiltered(&kept).pileup(contig, position).unwrap(),
            ));
        }
    }
    assert_eq!(swept, one_by_one);
    assert_eq!(swept, unfiltered);
    assert!(swept.iter().map(Vec::len).sum::<usize>() > 0);
}

#[test]
fn test_builder_reads_a_bam_stream() {
    let bytes = bam_bytes(HEADER, &[read("r", 10, "4M", "ACGT"), unmapped("u")]);
    let mut reader = bam::io::Reader::new(&bytes[..]);
    let header = reader.read_header().unwrap();
    let mut builder = StreamingPileupBuilder::new(reader, &header).unwrap();
    let mut columns = builder.columns("chr1", 9, 15).unwrap();
    let mut swept = Vec::new();
    while let Some(pileup) = columns.next_pileup() {
        swept.push(bases(&pileup.unwrap()));
    }
    assert_eq!(swept, ["", "A", "C", "G", "T", ""]);
}

#[test]
fn test_a_truncated_stream_is_an_error() {
    let bytes = bam_bytes(
        HEADER,
        &[
            read("r", 10, "4M", "ACGT"),
            read("s", 20, "4M", "ACGT"),
            read("t", 30, "4M", "ACGT"),
        ],
    );
    let mut decoded = Vec::new();
    io::Read::read_to_end(
        &mut noodles::bgzf::io::Reader::new(&bytes[..]),
        &mut decoded,
    )
    .unwrap();
    decoded.truncate(decoded.len() - 10);
    let mut reader = bam::io::Reader::from(&decoded[..]);
    let header = reader.read_header().unwrap();
    let mut builder = StreamingPileupBuilder::new(reader, &header).unwrap();
    assert_eq!(builder.pileup("chr1", 10).unwrap().len(), 1);
    assert!(matches!(builder.pileup("chr1", 20), Err(Error::Io(_))));
}

#[test]
fn test_the_kind_of_each_entry_is_named_as_in_python() {
    let names: Vec<&str> = [
        EntryKind::Base,
        EntryKind::Deletion,
        EntryKind::Insertion,
        EntryKind::Skip,
    ]
    .map(EntryKind::as_str)
    .into();
    assert_eq!(names, ["base", "deletion", "insertion", "skip"]);
}

#[test]
fn test_a_pileup_is_built_again_after_a_failed_advance() {
    let reads = [
        read("one", 10, "4M", "ACGT"),
        read("two", 20, "4M", "GGGG"),
        read("three", 5, "4M", "TTTT"),
    ];
    let mut builder = builder(&reads);
    assert_eq!(names(&builder.pileup("chr1", 10).unwrap()), ["one"]);
    for _ in 0..2 {
        assert!(matches!(
            builder.pileup("chr1", 21),
            Err(Error::OutOfOrder { name }) if name == "three"
        ));
    }

    let reads = [read("one", 10, "4M", "ACGT"), read("two", 20, "4M", "GGGG")];
    let tapped = log();
    let sink = Arc::clone(&tapped);
    let mut failed = false;
    let mut builder = super::builder(&reads).tap(move |record| {
        sink.lock().unwrap().push(name(&record));
        if failed {
            Ok(())
        } else {
            failed = true;
            Err(io::Error::other("disk full"))
        }
    });
    assert_eq!(names(&builder.pileup("chr1", 10).unwrap()), ["one"]);
    assert!(matches!(builder.pileup("chr1", 20), Err(Error::Io(_))));
    let retried = builder.pileup("chr1", 20).unwrap();
    assert_eq!(
        (names(&retried), bases(&retried)),
        (vec!["two".into()], "G".into())
    );
    builder.close().unwrap();
    assert_eq!(logged(&tapped), ["one", "two"]);
}

#[test]
fn test_a_read_that_fails_to_be_accepted_still_reaches_the_tap() {
    let (header, mut records) = records(
        HEADER,
        &[
            read("before", 10, "4M", "ACGT"),
            read("after", 30, "4M", "ACGT"),
        ],
    );
    records.insert(1, raw_record("bad", 11, &[4 << 4], 5));
    let evicted = log();
    let mut builder = StreamingPileupBuilder::new(source(records.clone()), &header)
        .unwrap()
        .tap(tap_into(&evicted));
    assert!(matches!(
        builder.pileup("chr1", 12),
        Err(Error::InvalidRecord { name, .. }) if name == "bad"
    ));
    assert_eq!(names(&builder.pileup("chr1", 31).unwrap()), ["after"]);
    builder.close().unwrap();
    assert_eq!(logged(&evicted), ["before", "bad", "after"]);

    let mut builder = StreamingPileupBuilder::new(source(records), &header).unwrap();
    assert!(builder.pileup("chr1", 12).is_err());
    let held: Vec<String> = builder
        .active
        .iter()
        .chain(builder.next.iter())
        .map(|&index| name(&builder.slots[index as usize].record))
        .collect();
    assert_eq!(held, ["before"]);
    assert_eq!(builder.free.len(), builder.slots.len() - 1);
}

#[test]
fn test_closing_again_after_a_tap_error_hands_over_the_rest() {
    let reads = [
        read("one", 10, "4M", "ACGT"),
        read("two", 20, "4M", "ACGT"),
        read("three", 30, "4M", "ACGT"),
    ];
    let tapped = log();
    let sink = Arc::clone(&tapped);
    let mut failed = false;
    let mut builder = builder(&reads).tap(move |record| {
        let record_name = name(&record);
        if record_name == "one" && !failed {
            failed = true;
            return Err(io::Error::other("disk full"));
        }
        sink.lock().unwrap().push(record_name);
        Ok(())
    });
    builder.pileup("chr1", 10).unwrap();
    assert!(matches!(builder.close(), Err(Error::Io(_))));
    assert!(matches!(builder.pileup("chr1", 40), Err(Error::Closed)));
    builder.close().unwrap();
    builder.close().unwrap();
    assert_eq!(logged(&tapped), ["two", "three"]);
}

#[test]
fn test_a_cigar_whose_query_offsets_would_overflow_is_refused() {
    let longest = 0x0FFF_FFFF_u32;
    let mut cigar = vec![(longest << 4) | 4; 15];
    cigar.extend([(longest << 4) | 1, longest << 4]);
    let (header, _) = records(HEADER, &[]);
    let source = source(vec![raw_record("r", 10, &cigar, 4)]);
    let mut builder = StreamingPileupBuilder::new(source, &header).unwrap();
    assert!(matches!(
        builder.pileup("chr1", 24),
        Err(Error::InvalidRecord { ref name, ref source })
            if name == "r" && source.to_string() == "CIGAR and query sequence lengths differ"
    ));
}
