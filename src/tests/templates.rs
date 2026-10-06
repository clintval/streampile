use super::{Read, read, unfiltered};
use crate::testing::{Pair, SamBuilder, Strand};
use crate::{AgreementStrategy, DisagreementStrategy, EntryKind, PileupTemplate};

const AGREEMENTS: [AgreementStrategy; 3] = [
    AgreementStrategy::Consensus,
    AgreementStrategy::MaxQual,
    AgreementStrategy::PassThrough,
];

const DISAGREEMENTS: [DisagreementStrategy; 3] = [
    DisagreementStrategy::Consensus,
    DisagreementStrategy::MaskBoth,
    DisagreementStrategy::MaskLowerQual,
];

/// A template's name, kind, base, and quality.
type Called = (String, EntryKind, Option<char>, Option<u8>);

fn called(template: &PileupTemplate<'_>) -> Called {
    (
        template.name().to_string(),
        template.kind(),
        template.base().map(char::from),
        template.quality(),
    )
}

/// The templates of the reads at a position, piled up without a filter or a quality floor.
fn templates_at(
    reads: &[Read],
    position: usize,
    agreement: AgreementStrategy,
    disagreement: DisagreementStrategy,
) -> Vec<Called> {
    let mut builder = unfiltered(reads).min_base_quality(0);
    let pileup = builder.pileup("chr1", position).unwrap();
    pileup
        .templates(agreement, disagreement)
        .iter()
        .map(called)
        .collect()
}

/// The base and quality of the one template at a position.
fn base_at(
    reads: &[Read],
    position: usize,
    agreement: AgreementStrategy,
    disagreement: DisagreementStrategy,
) -> (Option<char>, Option<u8>) {
    let templates = templates_at(reads, position, agreement, disagreement);
    assert_eq!(templates.len(), 1, "{reads:?}");
    (templates[0].2, templates[0].3)
}

/// Two overlapping mates of an FR pair whose bases at 0-based 12 are given, at these qualities.
fn mates(base1: char, quality1: u8, base2: char, quality2: u8) -> Vec<Read> {
    let bases = |base: char| format!("AA{base}A");
    vec![
        read("t", 10, "4M", &bases(base1))
            .flag(99)
            .quals(&[30, 30, quality1, 30]),
        read("t", 10, "4M", &bases(base2))
            .flag(147)
            .quals(&[30, 30, quality2, 30]),
    ]
}

#[test]
fn test_templates_group_reads_by_name_in_the_order_of_their_first_entries() {
    let reads = [
        read("q3", 50, "50M", &"A".repeat(50)).flag(99),
        read("q1", 100, "50M", &"C".repeat(50)).flag(99),
        read("q2", 100, "50M", &"G".repeat(50)).flag(147),
        read("q3", 100, "50M", &"T".repeat(50)).flag(147),
        read("q1", 110, "50M", &"C".repeat(50)).flag(147),
        read("q2", 110, "50M", &"G".repeat(50)).flag(99),
    ];
    let mut builder = super::builder(&reads);
    let pileup = builder.pileup("chr1", 125).unwrap();
    assert_eq!(pileup.unfiltered_depth(), 5);
    let templates = pileup.templates(
        AgreementStrategy::default(),
        DisagreementStrategy::default(),
    );
    let reads: Vec<usize> = templates
        .iter()
        .map(|template| template.entries().count())
        .collect();
    assert_eq!(reads, [2, 2, 1]);
    let called: Vec<Called> = templates.iter().map(called).collect();
    assert_eq!(
        called,
        [
            ("q1".into(), EntryKind::Base, Some('C'), Some(80)),
            ("q2".into(), EntryKind::Base, Some('G'), Some(80)),
            ("q3".into(), EntryKind::Base, Some('T'), Some(40)),
        ]
    );
}

#[test]
fn test_the_strategies_default_to_consensus_as_in_fgbio_and_fgumi() {
    assert_eq!(AgreementStrategy::default(), AgreementStrategy::Consensus);
    assert_eq!(
        DisagreementStrategy::default(),
        DisagreementStrategy::Consensus
    );
}

#[test]
fn test_agreeing_bases_make_the_quality_of_the_agreement_strategy() {
    let expected = [
        (AgreementStrategy::Consensus, (30, 35), 65),
        (AgreementStrategy::Consensus, (60, 50), 93),
        (AgreementStrategy::MaxQual, (30, 35), 35),
        (AgreementStrategy::MaxQual, (60, 50), 60),
        (AgreementStrategy::PassThrough, (30, 35), 35),
        (AgreementStrategy::PassThrough, (60, 50), 60),
    ];
    for (agreement, (quality1, quality2), quality) in expected {
        let reads = mates('C', quality1, 'C', quality2);
        for disagreement in DISAGREEMENTS {
            assert_eq!(
                base_at(&reads, 12, agreement, disagreement),
                (Some('C'), Some(quality)),
                "{agreement:?} of Q{quality1} and Q{quality2}"
            );
        }
    }
}

#[test]
fn test_disagreeing_bases_make_the_base_and_quality_of_the_disagreement_strategy() {
    let expected = [
        (
            DisagreementStrategy::Consensus,
            (30, 20),
            ('A', 10),
            ('C', 10),
        ),
        (
            DisagreementStrategy::Consensus,
            (21, 20),
            ('A', 2),
            ('C', 2),
        ),
        (
            DisagreementStrategy::Consensus,
            (30, 30),
            ('N', 2),
            ('N', 2),
        ),
        (DisagreementStrategy::MaskBoth, (30, 20), ('N', 2), ('N', 2)),
        (DisagreementStrategy::MaskBoth, (30, 30), ('N', 2), ('N', 2)),
        (
            DisagreementStrategy::MaskLowerQual,
            (30, 20),
            ('A', 30),
            ('C', 30),
        ),
        (
            DisagreementStrategy::MaskLowerQual,
            (1, 0),
            ('A', 1),
            ('C', 1),
        ),
        (
            DisagreementStrategy::MaskLowerQual,
            (30, 30),
            ('N', 2),
            ('N', 2),
        ),
    ];
    for (disagreement, (high, low), first_higher, second_higher) in expected {
        for agreement in AGREEMENTS {
            let first = base_at(&mates('A', high, 'C', low), 12, agreement, disagreement);
            let second = base_at(&mates('A', low, 'C', high), 12, agreement, disagreement);
            assert_eq!(
                (first, second),
                (
                    (Some(first_higher.0), Some(first_higher.1)),
                    (Some(second_higher.0), Some(second_higher.1))
                ),
                "{disagreement:?} of Q{high} and Q{low}"
            );
        }
    }
}

#[test]
fn test_a_no_call_leaves_the_other_reads_base_at_its_own_quality() {
    for agreement in AGREEMENTS {
        for disagreement in DISAGREEMENTS {
            let called = |reads: &[Read]| base_at(reads, 12, agreement, disagreement);
            assert_eq!(called(&mates('N', 40, 'A', 20)), (Some('A'), Some(20)));
            assert_eq!(called(&mates('G', 20, 'N', 40)), (Some('G'), Some(20)));
            assert_eq!(called(&mates('N', 10, 'N', 30)), (Some('N'), Some(30)));
        }
    }
}

#[test]
fn test_a_read_with_a_deletion_or_a_skip_holds_no_base() {
    let base = read("t", 10, "4M", "ACGT")
        .flag(99)
        .quals(&[30, 30, 15, 30]);
    let deletion = read("t", 10, "2M1D2M", "ACTT")
        .flag(147)
        .quals(&[30, 30, 25, 30]);
    let skip = read("t", 10, "2M1N2M", "ACTT").flag(147);
    let default = (
        AgreementStrategy::default(),
        DisagreementStrategy::default(),
    );
    let called = |reads: &[Read]| templates_at(reads, 12, default.0, default.1);
    let expected = |kind, base, quality| vec![("t".to_owned(), kind, base, quality)];
    assert_eq!(
        called(&[base.clone(), deletion.clone()]),
        expected(EntryKind::Base, Some('G'), Some(15))
    );
    assert_eq!(
        called(&[base, skip.clone()]),
        expected(EntryKind::Base, Some('G'), Some(15))
    );
    let deletion1 = read("t", 10, "2M1D2M", "ACTT")
        .flag(99)
        .quals(&[30, 30, 20, 30]);
    assert_eq!(
        called(&[deletion1.clone(), deletion]),
        expected(EntryKind::Deletion, None, Some(25))
    );
    assert_eq!(
        called(&[deletion1, skip.clone()]),
        expected(EntryKind::Deletion, None, Some(20))
    );
    let skip1 = read("t", 10, "2M1N2M", "ACTT").flag(99);
    assert_eq!(
        called(&[skip1, skip]),
        expected(EntryKind::Skip, None, None)
    );
}

#[test]
fn test_insertion_entries_are_no_part_of_a_template() {
    let reads = [
        read("t", 10, "3M2I3M", "ACGTTACG").flag(99),
        read("t", 10, "6M", "ACGACG").flag(147),
        read("opens", 13, "1I3M", "TACG"),
    ];
    let mut builder = unfiltered(&reads);
    let pileup = builder.pileup("chr1", 12).unwrap();
    assert_eq!(pileup.len(), 4);
    let templates = pileup.templates(
        AgreementStrategy::default(),
        DisagreementStrategy::default(),
    );
    assert_eq!(templates.len(), 1);
    let kinds: Vec<EntryKind> = templates[0].entries().map(|entry| entry.kind()).collect();
    assert_eq!(kinds, [EntryKind::Base, EntryKind::Base]);
    assert_eq!(
        called(&templates[0]),
        ("t".into(), EntryKind::Base, Some('G'), Some(80))
    );
}

/// The templates of the reads at a position under a quality floor and the default strategies.
fn floored(reads: &[Read], position: usize, floor: u8) -> Vec<Called> {
    let mut builder = unfiltered(reads).min_base_quality(floor);
    let pileup = builder.pileup("chr1", position).unwrap();
    pileup
        .templates(
            AgreementStrategy::default(),
            DisagreementStrategy::default(),
        )
        .iter()
        .map(called)
        .collect()
}

#[test]
fn test_a_read_under_the_floor_does_not_vote() {
    let reads = mates('C', 10, 'C', 10);
    let mut builder = unfiltered(&reads).min_base_quality(13);
    let pileup = builder.pileup("chr1", 12).unwrap();
    assert_eq!(pileup.filtered_depth(), 0);
    let templates = pileup.templates(
        AgreementStrategy::default(),
        DisagreementStrategy::default(),
    );
    assert_eq!(templates[0].entries().count(), 2);
    assert_eq!(
        called(&templates[0]),
        ("t".into(), EntryKind::Base, None, None)
    );
    assert!(!templates[0].passes(0));
    assert_eq!(
        floored(&reads, 12, 0),
        [("t".into(), EntryKind::Base, Some('C'), Some(20))]
    );
    assert_eq!(
        floored(&mates('C', 10, 'C', 30), 12, 13),
        [("t".into(), EntryKind::Base, Some('C'), Some(30))]
    );
}

#[test]
fn test_a_mate_under_the_floor_neither_masks_nor_lowers_the_other_mates_base() {
    let reads = mates('A', 35, 'C', 5);
    for disagreement in DISAGREEMENTS {
        let mut builder = unfiltered(&reads).min_base_quality(13);
        let pileup = builder.pileup("chr1", 12).unwrap();
        let templates = pileup.templates(AgreementStrategy::default(), disagreement);
        assert_eq!(
            (templates[0].base(), templates[0].quality()),
            (Some(b'A'), Some(35)),
            "{disagreement:?}"
        );
    }
    let masked = templates_at(
        &reads,
        12,
        AgreementStrategy::default(),
        DisagreementStrategy::MaskBoth,
    );
    assert_eq!(masked[0].2, Some('N'));
}

#[test]
fn test_a_deletion_votes_by_the_quality_of_its_next_base() {
    let base = read("t", 10, "4M", "ACGT").flag(99).quals(&[30, 30, 5, 30]);
    let deletion = read("t", 10, "2M1D2M", "ACTT")
        .flag(147)
        .quals(&[30, 30, 25, 30]);
    assert_eq!(
        floored(&[base.clone(), deletion], 12, 13),
        [("t".into(), EntryKind::Deletion, None, Some(25))]
    );
    let shallow = read("t", 10, "2M1D2M", "ACTT")
        .flag(147)
        .quals(&[30, 30, 5, 30]);
    assert_eq!(
        floored(&[base, shallow], 12, 13),
        [("t".into(), EntryKind::Base, None, None)]
    );
}

#[test]
fn test_more_than_two_reads_of_a_name_are_called_in_input_order() {
    let reads = [
        read("t", 10, "4M", "ACGT").flag(99).quals(&[20; 4]),
        read("t", 10, "4M", "ACGT").flag(2147).quals(&[20; 4]),
        read("t", 10, "4M", "ACAT").flag(147).quals(&[30; 4]),
    ];
    let called = |agreement, disagreement| base_at(&reads, 12, agreement, disagreement);
    assert_eq!(
        called(
            AgreementStrategy::Consensus,
            DisagreementStrategy::Consensus
        ),
        (Some('G'), Some(10))
    );
    assert_eq!(
        called(
            AgreementStrategy::MaxQual,
            DisagreementStrategy::MaskLowerQual
        ),
        (Some('A'), Some(30))
    );
    let mut builder = unfiltered(&reads);
    let pileup = builder.pileup("chr1", 12).unwrap();
    let templates = pileup.templates(
        AgreementStrategy::default(),
        DisagreementStrategy::default(),
    );
    assert_eq!(templates[0].entries().count(), 3);
}

#[test]
fn test_a_templates_strand_is_its_first_reads() {
    let forward = read("t", 10, "4M", "ACGT").flag(99);
    let reverse = read("t", 10, "4M", "ACGT").flag(147);
    let cases = [
        (vec![forward.clone(), reverse.clone()], false),
        (
            vec![reverse.clone().flag(83), forward.clone().flag(163)],
            true,
        ),
        (vec![reverse], false),
        (vec![forward.flag(163)], true),
        (vec![read("t", 10, "4M", "ACGT").flag(16)], true),
        (vec![read("t", 10, "4M", "ACGT")], false),
    ];
    for (reads, reverse) in cases {
        let mut builder = unfiltered(&reads);
        let pileup = builder.pileup("chr1", 12).unwrap();
        let templates = pileup.templates(
            AgreementStrategy::default(),
            DisagreementStrategy::default(),
        );
        assert_eq!(templates[0].is_reverse(), reverse, "{reads:?}");
    }
}

/// The templates' distances to both of their ends at each 0-based position of built pairs.
fn template_ends(pair: Pair, positions: &[usize]) -> Vec<(Option<usize>, Option<usize>)> {
    let mut builder = SamBuilder::new();
    builder.add_pair(pair);
    let mut pileups = builder.to_pileup_builder();
    positions
        .iter()
        .map(|&position| {
            let pileup = pileups.pileup("chr1", position).unwrap();
            let templates = pileup.templates(
                AgreementStrategy::default(),
                DisagreementStrategy::default(),
            );
            let template = &templates[0];
            (
                template.five_prime_distance().unwrap(),
                template.template_end_distance().unwrap(),
            )
        })
        .collect()
}

#[test]
fn test_a_templates_distances_are_its_first_reads_and_else_its_second_reads() {
    let pair = Pair {
        cigar1: Some("10M".into()),
        cigar2: Some("10M".into()),
        ..Pair::at(101, 106)
    };
    let expected = [(0, 14), (4, 10), (5, 9), (9, 5), (10, 4), (14, 0)];
    assert_eq!(
        template_ends(pair, &[100, 104, 105, 109, 110, 114]),
        expected.map(|(five, end)| (Some(five), Some(end)))
    );
    let deleted = Pair {
        cigar1: Some("2M1D7M".into()),
        cigar2: Some("10M".into()),
        ..Pair::at(101, 101)
    };
    assert_eq!(template_ends(deleted, &[102]), [(Some(2), Some(7))]);
}

#[test]
fn test_a_fragment_has_no_template_end() {
    let reads = [read("f", 10, "4M", "ACGT").flag(16)];
    let mut builder = unfiltered(&reads);
    let pileup = builder.pileup("chr1", 11).unwrap();
    let templates = pileup.templates(
        AgreementStrategy::default(),
        DisagreementStrategy::default(),
    );
    let template = &templates[0];
    assert_eq!(
        (
            template.five_prime_distance().unwrap(),
            template.template_end_distance().unwrap(),
            template.is_no_call(),
            template.passes(40),
        ),
        (Some(2), None, false, true)
    );
}

#[test]
fn test_a_pileup_without_a_base_deletion_or_skip_has_no_templates() {
    let reads = [read("opens", 13, "1I3M", "TACG")];
    let mut builder = unfiltered(&reads);
    let empty = builder.pileup("chr1", 5).unwrap();
    assert!(
        empty
            .templates(
                AgreementStrategy::default(),
                DisagreementStrategy::default()
            )
            .is_empty()
    );
    let insertion = builder.pileup("chr1", 12).unwrap();
    assert_eq!(insertion.len(), 1);
    assert!(
        insertion
            .templates(
                AgreementStrategy::default(),
                DisagreementStrategy::default()
            )
            .is_empty()
    );
}

#[test]
fn test_a_template_of_a_pair_that_is_not_fr_has_no_template_end() {
    let outward = Pair {
        cigar1: Some("10M".into()),
        cigar2: Some("10M".into()),
        strand1: Strand::Minus,
        strand2: Strand::Plus,
        ..Pair::at(101, 201)
    };
    assert_eq!(
        template_ends(outward, &[105, 205]),
        [(Some(4), None), (None, None)]
    );
}
