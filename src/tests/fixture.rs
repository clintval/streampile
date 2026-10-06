use std::fs::File;

use noodles::sam::alignment::record::Flags;
use noodles::{bam, bgzf};

use super::{bases, names};
use crate::{AgreementStrategy, AuxValue, DisagreementStrategy, StreamingPileupBuilder};

const READS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/reads.bam");

type Reader = bam::io::Reader<bgzf::io::Reader<File>>;

fn open<'f>() -> StreamingPileupBuilder<'f, Reader> {
    let mut reader = bam::io::reader::Builder
        .build_from_path(READS)
        .expect("the fixture BAM");
    let header = reader.read_header().expect("the fixture header");
    StreamingPileupBuilder::new(reader, &header).expect("a coordinate-sorted fixture")
}

#[test]
fn test_the_fixture_piles_up_as_the_readme_shows() {
    let mut builder = open().min_base_quality(30);
    let first = builder.pileup("chr1", 10).unwrap();
    assert_eq!(
        (first.filtered_depth(), bases(&first)),
        (4, "ATGA".to_owned())
    );
    let second = builder.pileup("chr1", 12).unwrap();
    assert_eq!(
        (second.filtered_depth(), bases(&second)),
        (3, "GGG".to_owned())
    );
}

#[test]
fn test_the_fixture_pair_is_one_template_whose_mates_agree_and_then_disagree() {
    let paired = |record: &bam::Record| record.flags().is_segmented();
    let mut builder = open().read_filter(paired);
    let mut called = |position, agreement, disagreement| {
        let pileup = builder.pileup("chr1", position).unwrap();
        let templates = pileup.templates(agreement, disagreement);
        assert_eq!(templates.len(), 1);
        let template = &templates[0];
        assert_eq!(template.entries().count(), 2);
        (template.base().map(char::from), template.quality())
    };
    assert_eq!(
        called(
            50,
            AgreementStrategy::Consensus,
            DisagreementStrategy::Consensus
        ),
        (Some('T'), Some(80))
    );
    let disagreeing: Vec<(Option<char>, Option<u8>)> = [
        DisagreementStrategy::Consensus,
        DisagreementStrategy::MaskBoth,
        DisagreementStrategy::MaskLowerQual,
    ]
    .into_iter()
    .map(|strategy| called(52, AgreementStrategy::Consensus, strategy))
    .collect();
    assert_eq!(
        disagreeing,
        [
            (Some('G'), Some(20)),
            (Some('N'), Some(2)),
            (Some('G'), Some(40))
        ]
    );
}

#[test]
fn test_the_fixture_sweeps_a_span() {
    let mut builder = open();
    let mut columns = builder.columns("chr1", 0, 20).unwrap();
    let mut depth = 0;
    while let Some(pileup) = columns.next_pileup() {
        depth += pileup.unwrap().unfiltered_depth();
    }
    assert_eq!(depth, 73);
}

#[test]
fn test_the_fixture_reads_aux_fields_in_place() {
    let mut builder = open()
        .exclude_flags(Flags::empty())
        .index_aux_tags([*b"RG"]);
    let pileup = builder.pileup("chr1", 56).unwrap();
    assert_eq!(names(&pileup), ["pair"]);
    let entry = pileup.get(0).unwrap();
    assert_eq!(entry.aux(*b"RG").unwrap(), Some(AuxValue::String(b"A")));
    assert_eq!(entry.aux(*b"NM").unwrap(), None);
    let mut builder = open();
    let pileup = builder.pileup("chr1", 56).unwrap();
    assert_eq!(
        pileup.get(0).unwrap().aux(*b"RG").unwrap(),
        Some(AuxValue::String(b"A"))
    );
}

#[test]
fn test_the_fixture_pair_measures_both_fragment_ends_from_its_mate_cigars() {
    let mut builder = open();
    let pileup = builder.pileup("chr1", 45).unwrap();
    let distances: Vec<(u16, Option<usize>, Option<usize>)> = pileup
        .iter()
        .map(|entry| {
            (
                entry.flags().bits(),
                entry.five_prime_distance(),
                entry.template_end_distance().unwrap(),
            )
        })
        .collect();
    assert_eq!(
        distances,
        [
            (0, Some(5), None),
            (99, Some(5), Some(14)),
            (147, Some(14), Some(5))
        ]
    );
}
