//! Tests of the FR rule of htsjdk 5.0.0's `SamPairUtil.getPairOrientation` (samtools/htsjdk#1771),
//! ported from htsjdk's `SamPairUtilTest`, fulcrumgenomics/fgumi#1022, and fgbio's
//! `SamRecordTest` and `CodecConsensusCallerTest`, each cited where it is ported.
//!
//! Positions in the pairs built here are 1-based, as in htsjdk and fgbio, and positions given to
//! `template_end_distance` are 0-based. The cases of fgumi#1022 that test the forward record's
//! TLEN arm are not ported, since the mate's end is read from `MC` and TLEN is never read.

use noodles::sam;
use noodles::sam::alignment::RecordBuf;

use crate::testing::{Frag, Pair, SamBuilder, Strand};
use crate::{Error, is_fr_pair, template_end_distance};

/// A pair of reads of all `A`s with these CIGARs at these 1-based starts and strands, as fgbio's
/// `SamBuilder.addPair` and htsjdk's `SamPairUtil.setMateInfo` make them.
fn pair(
    start1: usize,
    cigar1: &str,
    reverse1: bool,
    start2: usize,
    cigar2: &str,
    reverse2: bool,
) -> (sam::Header, Vec<RecordBuf>) {
    let length = |cigar: &str| {
        noodles::sam::record::Cigar::new(cigar.as_bytes())
            .iter()
            .map(Result::unwrap)
            .filter(|op| op.kind().consumes_read())
            .map(noodles::sam::alignment::record::cigar::Op::len)
            .sum::<usize>()
    };
    let strand = |reverse: bool| if reverse { Strand::Minus } else { Strand::Plus };
    let mut builder = SamBuilder::new();
    let records = builder.add_pair(Pair {
        name: Some("q".into()),
        bases1: Some("A".repeat(length(cigar1))),
        bases2: Some("A".repeat(length(cigar2))),
        cigar1: Some(cigar1.into()),
        cigar2: Some(cigar2.into()),
        strand1: strand(reverse1),
        strand2: strand(reverse2),
        ..Pair::at(start1, start2)
    });
    (builder.header().clone(), records)
}

/// Whether each read of a pair is a read of an FR pair.
fn fr(header: &sam::Header, records: &[RecordBuf]) -> Vec<bool> {
    records
        .iter()
        .map(|record| is_fr_pair(record, header).unwrap())
        .collect()
}

/// Each read's template-end distance at a 0-based position.
fn ends(header: &sam::Header, records: &[RecordBuf], position: usize) -> Vec<Option<usize>> {
    records
        .iter()
        .map(|record| template_end_distance(record, header, position).unwrap())
        .collect()
}

/// A test's name, each read's 1-based start, length, and whether it is reverse, and whether the
/// pair is FR.
type Vector = (&'static str, usize, usize, bool, usize, usize, bool, bool);

/// htsjdk 5.0.0 `SamPairUtilTest`'s orientation vectors.
#[rustfmt::skip]
const VECTORS: [Vector; 16] = [
    ("normal innie", 1, 100, false, 500, 100, true, true),
    ("overlapping innie", 1, 100, false, 50, 100, true, true),
    ("second end enclosed innie", 1, 100, false, 50, 50, true, true),
    ("first end enclosed innie", 1, 50, false, 1, 100, true, true),
    ("completely overlapping innie", 1, 100, false, 1, 100, true, true),
    ("normal outie", 1, 100, true, 500, 100, false, false),
    ("nojump outie", 1, 100, true, 101, 100, false, false),
    ("forward tandem", 1, 100, true, 500, 100, true, false),
    ("reverse tandem", 1, 100, false, 500, 100, false, false),
    ("overlapping forward tandem", 1, 100, true, 50, 100, true, false),
    ("overlapping reverse tandem", 1, 100, false, 50, 100, false, false),
    ("second end enclosed forward tandem", 1, 100, true, 50, 50, true, false),
    ("second end enclosed reverse tandem", 1, 100, false, 50, 50, false, false),
    ("first end enclosed forward tandem", 1, 50, true, 1, 100, true, false),
    ("first end enclosed reverse tandem", 1, 50, false, 1, 100, false, false),
    ("dovetail 5' tie", 100, 100, false, 1, 100, true, true),
];

/// htsjdk 5.0.0 `SamPairUtilTest.testGetPairOrientation`: its data provider's vectors, each
/// asserted on both ends after `setMateInfo`. FR is `true`, and RF and tandem are `false`.
#[test]
fn test_htsjdk_pair_orientation_vectors() {
    for (name, start1, length1, reverse1, start2, length2, reverse2, expected) in VECTORS {
        let (header, records) = pair(
            start1,
            &format!("{length1}M"),
            reverse1,
            start2,
            &format!("{length2}M"),
            reverse2,
        );
        assert_eq!(fr(&header, &records), [expected, expected], "{name}");
    }
}

/// htsjdk 5.0.0 `SamPairUtilTest` "dovetail 5' tie": the forward read `100M` at 100 and the
/// reverse read `100M` at 1 share only position 100, where both 5′ ends lie, so each read is 0
/// bases from its mate's 5′ end there.
#[test]
fn test_htsjdk_dovetail_five_prime_tie_is_fr_from_both_reads() {
    let (header, records) = pair(100, "100M", false, 1, "100M", true);
    assert_eq!(fr(&header, &records), [true, true]);
    assert_eq!(ends(&header, &records, 99), [Some(0), Some(0)]);
}

/// htsjdk 5.0.0 `SamPairUtilTest.testGetPairOrientationSymmetryForSoftClippedDovetail`: R1
/// forward `90S10M` at 100 (5′ at 100, aligned end 109) and R2 reverse `10M90S` at 91 (aligned
/// end and 5′ at 100) are FR from both reads. At 100, the one position both reads cover, each
/// read's template end counts its mate's 90 soft-clipped bases.
#[test]
fn test_htsjdk_soft_clipped_dovetail_is_fr_from_both_reads() {
    let (header, records) = pair(100, "90S10M", false, 91, "10M90S", true);
    assert_eq!(fr(&header, &records), [true, true]);
    assert_eq!(ends(&header, &records, 99), [Some(90), Some(90)]);
}

/// fgumi#1022 `test_is_primary_fr_pair_raw_coincident_five_prime_dovetail`, on fgbio's
/// `CodecConsensusCallerTest` pair "classify a dovetail FR pair as FR regardless of argument
/// order": forward `68S53M8S` at 96 (aligned 96 to 148) and reverse `28S48M53S` at 49 (aligned
/// 49 to 96), whose aligned 5′ ends meet at 96. Both reads, in either order of the pair, are FR,
/// and at 96 each read's template end counts its mate's soft-clipped bases past it, 53 and 68.
#[test]
fn test_fgumi_coincident_five_prime_dovetail_is_fr_in_either_order() {
    let (header, records) = pair(96, "68S53M8S", false, 49, "28S48M53S", true);
    assert_eq!(fr(&header, &records), [true, true]);
    let (header, swapped) = pair(49, "28S48M53S", true, 96, "68S53M8S", false);
    assert_eq!(fr(&header, &swapped), [true, true]);
    assert_eq!(ends(&header, &records, 95), [Some(53), Some(68)]);
}

/// fgumi#1022 `test_is_fr_pair_raw_strict_rf_stays_not_fr` and
/// `test_get_pair_orientation_one_past_five_prime_tie_is_rf`: one base past the tie, a forward
/// read starting right after its reverse mate's end, is RF from both reads.
#[test]
fn test_fgumi_one_past_the_five_prime_tie_is_rf() {
    let (header, records) = pair(100, "100M", true, 200, "100M", false);
    assert_eq!(fr(&header, &records), [false, false]);
    let (header, records) = pair(101, "100M", false, 1, "100M", true);
    assert_eq!(fr(&header, &records), [false, false]);
}

/// fgumi#1022 `test_zero_reference_span_reverse_read_at_mate_start_is_rf`: a reverse read with
/// no reference-consuming operator, `100S`, at its forward mate's start ends at 99, before its
/// start, so the pair is RF from both reads, and neither read has a template end.
#[test]
fn test_fgumi_zero_reference_span_reverse_read_at_its_mates_start_is_rf() {
    let (header, records) = pair(100, "100M", false, 100, "100S", true);
    assert_eq!(fr(&header, &records), [false, false]);
    assert_eq!(ends(&header, &records, 99), [None, None]);
}

/// fgumi#1022 widens the comparison to i64 so that extreme coordinates cannot overflow; here the
/// extremes are starts near the largest BAM position and spans past it, as TLEN is never read.
#[test]
fn test_fgumi_extreme_coordinates_do_not_overflow() {
    let long = "1M268435455N1M";
    let (header, records) = pair(2_100_000_000, long, false, 2_000_000_000, long, true);
    assert_eq!(fr(&header, &records), [true, true]);
    let (header, records) = pair(2_147_483_647, "1M", false, 2_147_483_000, long, true);
    assert_eq!(fr(&header, &records), [true, true]);
}

/// The template length is never read: an FR pair at a 5′ tie with bwa's TLEN of 0, or any other,
/// is FR from both reads, where htsjdk's TLEN fallback would call the forward read RF.
#[test]
fn test_the_template_length_is_never_read() {
    let (header, mut records) = pair(100, "100M", false, 1, "100M", true);
    for tlen in [0, -7, i32::MIN] {
        for record in &mut records {
            *record.template_length_mut() = tlen;
        }
        assert_eq!(fr(&header, &records), [true, true], "{tlen}");
    }
}

/// fgbio `SamRecordTest` "SamRecord.isFrPair should return false for reads that are not part of
/// a mapped pair".
#[test]
fn test_fgbio_reads_not_of_a_mapped_pair_are_not_fr() {
    let mut builder = SamBuilder::new().read_length(10).base_quality(20);
    let mut records = builder.add_frag(Frag::at(100));
    records.extend(builder.add_pair(Pair {
        unmapped2: true,
        ..Pair::at(100, 100)
    }));
    assert_eq!(fr(builder.header(), &records), [false, false, false]);
}

/// fgbio `SamRecordTest` "SamRecord.isFrPair should return false for pairs that are mapped to
/// different chromosomes".
#[test]
fn test_fgbio_pairs_on_different_contigs_are_not_fr() {
    let mut builder = SamBuilder::new().read_length(10).base_quality(20);
    let records: Vec<RecordBuf> = builder
        .add_pair(Pair::at(100, 200))
        .into_iter()
        .map(|record| SamBuilder::with_mate_reference_sequence_id(record, 1))
        .collect();
    assert_eq!(fr(builder.header(), &records), [false, false]);
}

/// fgbio `SamRecordTest` "SamRecord.isFrPair should return false for mapped FF, RR and RF pairs"
/// and "should return true on an actual FR pair".
#[test]
fn test_fgbio_only_an_actual_fr_pair_is_fr() {
    let strands = [
        (Strand::Plus, Strand::Plus, false),
        (Strand::Minus, Strand::Minus, false),
        (Strand::Minus, Strand::Plus, false),
        (Strand::Plus, Strand::Minus, true),
    ];
    for (strand1, strand2, expected) in strands {
        let mut builder = SamBuilder::new().read_length(10).base_quality(20);
        let records = builder.add_pair(Pair {
            strand1,
            strand2,
            ..Pair::at(100, 200)
        });
        assert_eq!(fr(builder.header(), &records), [expected, expected]);
    }
}

/// An outward-facing pair is not FR, so neither read has a template end, while a read of an FR
/// pair past its mate's 5′ end has none there and is still a read of an FR pair.
#[test]
fn test_not_fr_and_past_the_mates_end_are_told_apart_by_is_fr_pair() {
    let (header, records) = pair(100, "50M", true, 200, "50M", false);
    assert_eq!(fr(&header, &records), [false, false]);
    for position in [99, 120, 148, 199, 220, 248] {
        assert_eq!(
            ends(&header, &records, position),
            [None, None],
            "{position}"
        );
    }
    let (header, records) = pair(100, "100M", false, 60, "90M10S", true);
    assert_eq!(fr(&header, &records), [true, true]);
    let forward = &records[0];
    let ends: Vec<Option<usize>> = [157, 158, 159]
        .into_iter()
        .map(|position| template_end_distance(forward, &header, position).unwrap())
        .collect();
    assert_eq!(ends, [Some(1), Some(0), None]);
}

/// The forward read of an FR pair needs `MC` to find its mate's end, and the reverse read does
/// not; a pair that is not FR by its flags needs none.
#[test]
fn test_only_a_forward_read_of_an_fr_pair_needs_a_mate_cigar() {
    let (header, records) = pair(100, "50M", false, 150, "50M", true);
    let bare: Vec<RecordBuf> = records
        .into_iter()
        .map(SamBuilder::without_mate_cigar)
        .collect();
    let refused = is_fr_pair(&bare[0], &header);
    assert!(
        matches!(&refused, Err(Error::MissingMateCigar { name }) if name == "q"),
        "{refused:?}"
    );
    assert!(is_fr_pair(&bare[1], &header).unwrap());
    let (header, records) = pair(100, "50M", false, 150, "50M", false);
    let tandem = SamBuilder::without_mate_cigar(records[0].clone());
    assert!(!is_fr_pair(&tandem, &header).unwrap());
}
