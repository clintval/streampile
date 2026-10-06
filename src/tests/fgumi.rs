//! The template view against fgumi's `OverlappingBasesConsensusCaller`, which rewrites the bases
//! and qualities of a pair's overlapping mates as fgbio's `CallOverlappingConsensusBases` does.
//!
//! For random overlapping pairs, with mismatches, no-calls, qualities from 0 up, soft clips, and
//! indels in and out of the overlap, and for every pair of strategies, fgumi rewrites the pair and
//! the rewritten reads are piled up. At every position, the template that `templates()` calls from
//! the original pair must hold the base and quality of the rewritten reads' template: a called
//! base over a no-call, then the higher quality, then the first read, since the strategies leave
//! a read's base or quality unchanged only where the other read's is a no-call, is masked, or is
//! of lower or equal quality.
//!
//! The two differ by design where fgumi rewrites reads before any filtering, so pairs are piled up
//! with no filter and no quality floor, both reads mapped to one contig with stored qualities.

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::io;

use fgumi_consensus::{
    AgreementStrategy as FgumiAgreement, DisagreementStrategy as FgumiDisagreement,
    OverlappingBasesConsensusCaller,
};
use noodles::bam;
use noodles::sam;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::io::Write as _;
use noodles::sam::alignment::record::Flags;

use crate::testing::{Pair, SamBuilder, Strand};
use crate::{AgreementStrategy, DisagreementStrategy, EntryKind, Records, StreamingPileupBuilder};

const PAIRS: usize = 400;

const AGREEMENTS: [(AgreementStrategy, FgumiAgreement); 3] = [
    (AgreementStrategy::Consensus, FgumiAgreement::Consensus),
    (AgreementStrategy::MaxQual, FgumiAgreement::MaxQual),
    (AgreementStrategy::PassThrough, FgumiAgreement::PassThrough),
];

const DISAGREEMENTS: [(DisagreementStrategy, FgumiDisagreement); 3] = [
    (
        DisagreementStrategy::Consensus,
        FgumiDisagreement::Consensus,
    ),
    (DisagreementStrategy::MaskBoth, FgumiDisagreement::MaskBoth),
    (
        DisagreementStrategy::MaskLowerQual,
        FgumiDisagreement::MaskLowerQual,
    ),
];

const QUALITIES: [u8; 12] = [0, 1, 2, 3, 10, 20, 20, 30, 30, 37, 40, 45];

/// A xorshift generator, so that every run tests the same pairs.
struct Random(u64);

impl Random {
    fn below(&mut self, bound: usize) -> usize {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 % bound as u64) as usize
    }

    fn chance(&mut self, percent: usize) -> bool {
        self.below(100) < percent
    }

    fn base(&mut self) -> u8 {
        if self.chance(5) {
            b'N'
        } else {
            b"ACGT"[self.below(4)]
        }
    }

    fn quality(&mut self) -> u8 {
        QUALITIES[self.below(QUALITIES.len())]
    }
}

/// A read's CIGAR, bases, and qualities, with `aligned` bases from a 0-based start of a
/// reference, mostly its bases and sometimes an error or a no-call.
fn read(
    random: &mut Random,
    reference: &[u8],
    start: usize,
    aligned: usize,
) -> (String, String, Vec<u8>) {
    let (mut cigar, mut bases, mut qualities) = (String::new(), Vec::new(), Vec::new());
    let clip =
        |random: &mut Random, cigar: &mut String, bases: &mut Vec<u8>, qualities: &mut Vec<u8>| {
            if random.chance(40) {
                let length = 1 + random.below(5);
                write!(cigar, "{length}S").unwrap();
                for _ in 0..length {
                    bases.push(random.base());
                    qualities.push(random.quality());
                }
            }
        };
    clip(random, &mut cigar, &mut bases, &mut qualities);
    let (mut position, mut matched, mut done) = (start, 0, 0);
    while done < aligned {
        let roll = random.below(100);
        if matched > 0 && roll < 8 {
            write!(cigar, "{matched}M").unwrap();
            matched = 0;
            let length = 1 + random.below(3);
            if roll < 4 {
                write!(cigar, "{length}I").unwrap();
                for _ in 0..length {
                    bases.push(random.base());
                    qualities.push(random.quality());
                }
            } else {
                write!(cigar, "{length}D").unwrap();
                position += length;
            }
        } else {
            let base = if random.chance(7) {
                random.base()
            } else {
                reference[position]
            };
            bases.push(base);
            qualities.push(random.quality());
            position += 1;
            matched += 1;
            done += 1;
        }
    }
    write!(cigar, "{matched}M").unwrap();
    clip(random, &mut cigar, &mut bases, &mut qualities);
    let bases = String::from_utf8(bases).expect("ASCII bases");
    (cigar, bases, qualities)
}

/// A record as the raw BAM bytes fgumi rewrites, without the block size.
fn encode(header: &sam::Header, record: &RecordBuf) -> Vec<u8> {
    let mut writer = bam::io::Writer::from(Vec::new());
    writer.write_alignment_record(header, record).unwrap();
    writer.into_inner()[4..].to_vec()
}

fn decode(raw: &[u8]) -> bam::Record {
    let mut bytes = u32::try_from(raw.len()).unwrap().to_le_bytes().to_vec();
    bytes.extend_from_slice(raw);
    let mut record = bam::Record::default();
    bam::io::Reader::from(&bytes[..])
        .read_record(&mut record)
        .unwrap();
    record
}

/// The template of the bases of a pair's rewritten reads at a position.
fn rewritten_template(bases: &[(u8, u8)]) -> Option<(u8, u8)> {
    let called: Vec<(u8, u8)> = bases
        .iter()
        .copied()
        .filter(|(base, _)| *base != b'N')
        .collect();
    let candidates = if called.is_empty() { bases } else { &called };
    candidates
        .iter()
        .copied()
        .fold(None, |best: Option<(u8, u8)>, next| match best {
            Some(best) if best.1 >= next.1 => Some(best),
            _ => Some(next),
        })
}

/// What two bases of a pair at one position are, for counting what the test covers.
fn case(bases: &[(u8, u8)]) -> &'static str {
    match bases {
        [(a, _), (b, _)] if *a == b'N' || *b == b'N' => "a no-call",
        [(a, _), (b, _)] if a == b => "agreeing bases",
        [(_, qa), (_, qb)] if qa == qb => "disagreeing bases of equal quality",
        [_, _] => "disagreeing bases",
        _ => "one base",
    }
}

#[test]
fn test_templates_call_overlapping_mates_as_fgumi_does() {
    let mut random = Random(0x9E37_79B9_7F4A_7C15);
    let mut covered: BTreeMap<&'static str, usize> = BTreeMap::new();
    for index in 0..PAIRS {
        let reference: Vec<u8> = (0..600).map(|_| b"ACGT"[random.below(4)]).collect();
        let (start1, length1) = (100 + random.below(20), 20 + random.below(60));
        let (start2, length2) = (start1 + random.below(length1 + 10), 20 + random.below(60));
        let (cigar1, bases1, quals1) = read(&mut random, &reference, start1, length1);
        let (cigar2, bases2, quals2) = read(&mut random, &reference, start2, length2);
        let flipped = random.chance(30);
        let mut builder = SamBuilder::new();
        let records = builder.add_pair(Pair {
            name: Some(format!("p{index}")),
            bases1: Some(bases1),
            bases2: Some(bases2),
            quals1: Some(quals1),
            quals2: Some(quals2),
            cigar1: Some(cigar1.clone()),
            cigar2: Some(cigar2.clone()),
            strand1: if flipped { Strand::Minus } else { Strand::Plus },
            strand2: if flipped { Strand::Plus } else { Strand::Minus },
            ..Pair::at(start1 + 1, start2 + 1)
        });
        let header = builder.header().clone();
        for (agreement, fgumi_agreement) in AGREEMENTS {
            for (disagreement, fgumi_disagreement) in DISAGREEMENTS {
                let (mut raw1, mut raw2) =
                    (encode(&header, &records[0]), encode(&header, &records[1]));
                OverlappingBasesConsensusCaller::new(fgumi_agreement, fgumi_disagreement)
                    .call(&mut raw1, &mut raw2)
                    .unwrap();
                let mut rewritten = vec![decode(&raw1), decode(&raw2)];
                rewritten.sort_by_key(|record| record.alignment_start().unwrap().unwrap());
                let source = Records::new(rewritten.into_iter().map(Ok::<_, io::Error>));
                let mut theirs = StreamingPileupBuilder::new(source, &header)
                    .unwrap()
                    .exclude_flags(Flags::empty())
                    .min_base_quality(0);
                let mut ours = builder
                    .to_pileup_builder()
                    .exclude_flags(Flags::empty())
                    .min_base_quality(0);
                for position in start1.saturating_sub(10)..start2 + length2 + 30 {
                    let pileup = ours.pileup("chr1", position).unwrap();
                    let templates = pileup.templates(agreement, disagreement);
                    assert!(templates.len() <= 1);
                    let rewritten = theirs.pileup("chr1", position).unwrap();
                    let bases: Vec<(u8, u8)> = rewritten
                        .iter()
                        .filter(|entry| entry.kind() == EntryKind::Base)
                        .map(|entry| (entry.base().unwrap(), entry.quality().unwrap()))
                        .collect();
                    let expected = rewritten_template(&bases);
                    match templates.first() {
                        Some(template) if template.kind() == EntryKind::Base => {
                            let original: Vec<(u8, u8)> = template
                                .entries()
                                .filter(|entry| entry.kind() == EntryKind::Base)
                                .map(|entry| (entry.base().unwrap(), entry.quality().unwrap()))
                                .collect();
                            *covered.entry(case(&original)).or_default() += 1;
                            assert_eq!(
                                template.base().zip(template.quality()),
                                expected,
                                "{agreement:?} and {disagreement:?} at {position} of {cigar1} at \
                                 {start1} and {cigar2} at {start2}"
                            );
                        }
                        _ => assert_eq!(expected, None, "{position}"),
                    }
                }
            }
        }
    }
    for (case, least) in [
        ("agreeing bases", 5_000),
        ("disagreeing bases", 200),
        ("disagreeing bases of equal quality", 20),
        ("a no-call", 200),
        ("one base", 50_000),
    ] {
        let seen = covered.get(case).copied().unwrap_or(0);
        assert!(seen >= least, "{case}: {seen} of {covered:?}");
    }
}
