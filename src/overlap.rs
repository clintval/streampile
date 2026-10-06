//! How the reads of one template at one position become one base, with the strategies of fgbio's
//! `CallOverlappingConsensusBases` as fgumi implements them.

use std::cmp::Ordering;

use crate::pileup::EntryKind;

/// The highest quality the sum of two agreeing bases gets, as in fgumi.
const MAX_SUMMED_QUALITY: u16 = 93;

/// A no-call, the base of a masked base.
const NO_CALL: u8 = b'N';

/// The quality of a masked base.
const NO_CALL_QUALITY: u8 = 2;

/// How the quality of a template is made from two of its reads holding the same base.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum AgreementStrategy {
    /// The sum of the two qualities, at most 93.
    #[default]
    Consensus,
    /// The higher of the two qualities.
    MaxQual,
    /// Each read keeps its quality, so the template has the higher of the two.
    PassThrough,
}

impl AgreementStrategy {
    /// The strategy's name as the Python package spells it.
    pub fn as_str(self) -> &'static str {
        match self {
            AgreementStrategy::Consensus => "consensus",
            AgreementStrategy::MaxQual => "max_qual",
            AgreementStrategy::PassThrough => "pass_through",
        }
    }
}

/// How a template's base and quality are made from two of its reads holding different bases.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum DisagreementStrategy {
    /// The base of the higher quality, at the higher quality less the lower and at least 2, or an
    /// `N` at quality 2 when the qualities are equal.
    #[default]
    Consensus,
    /// An `N` at quality 2.
    MaskBoth,
    /// The base of the higher quality at its own quality, the other masked, or an `N` at quality 2
    /// when the qualities are equal.
    MaskLowerQual,
}

impl DisagreementStrategy {
    /// The strategy's name as the Python package spells it.
    pub fn as_str(self) -> &'static str {
        match self {
            DisagreementStrategy::Consensus => "consensus",
            DisagreementStrategy::MaskBoth => "mask_both",
            DisagreementStrategy::MaskLowerQual => "mask_lower_qual",
        }
    }
}

/// What one read holds at a position, or what a template's reads make of it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Observation {
    pub kind: EntryKind,
    pub base: Option<u8>,
    pub quality: Option<u8>,
}

/// What a template holds at a position, from what its reads hold there, in input order.
///
/// A read's base and the template's are called into one by the strategies, in input order. A
/// no-call is left alone, so the other read's base stands, and two no-calls are an `N` at the
/// higher quality. A read with a deletion or a skip holds no base, so a template with a base has it
/// at its own quality. A template with no base holds a deletion at the higher of the deletions'
/// qualities, or else a skip.
pub(crate) fn observe(
    reads: impl IntoIterator<Item = Observation>,
    agreement: AgreementStrategy,
    disagreement: DisagreementStrategy,
) -> Observation {
    let mut kind = EntryKind::Skip;
    let mut called: Option<(u8, u8)> = None;
    let mut deleted: Option<u8> = None;
    for read in reads {
        match read.kind {
            EntryKind::Base => {
                kind = EntryKind::Base;
                if let (Some(base), Some(quality)) = (read.base, read.quality) {
                    called = Some(match called {
                        Some(template) => call(template, (base, quality), agreement, disagreement),
                        None => (base, quality),
                    });
                }
            }
            EntryKind::Deletion => {
                if kind == EntryKind::Skip {
                    kind = EntryKind::Deletion;
                }
                deleted = deleted.max(read.quality);
            }
            EntryKind::Skip | EntryKind::Insertion => {}
        }
    }
    match kind {
        EntryKind::Base => Observation {
            kind,
            base: called.map(|(base, _)| base),
            quality: called.map(|(_, quality)| quality),
        },
        EntryKind::Deletion => Observation {
            kind,
            base: None,
            quality: deleted,
        },
        EntryKind::Skip | EntryKind::Insertion => Observation {
            kind: EntryKind::Skip,
            base: None,
            quality: None,
        },
    }
}

/// Two bases and their qualities called into one.
fn call(
    first: (u8, u8),
    second: (u8, u8),
    agreement: AgreementStrategy,
    disagreement: DisagreementStrategy,
) -> (u8, u8) {
    let ((base1, quality1), (base2, quality2)) = (first, second);
    match (base1 == NO_CALL, base2 == NO_CALL) {
        (true, true) => (NO_CALL, quality1.max(quality2)),
        (true, false) => second,
        (false, true) => first,
        (false, false) if base1 == base2 => {
            let quality = match agreement {
                AgreementStrategy::Consensus => {
                    (u16::from(quality1) + u16::from(quality2)).min(MAX_SUMMED_QUALITY) as u8
                }
                AgreementStrategy::MaxQual | AgreementStrategy::PassThrough => {
                    quality1.max(quality2)
                }
            };
            (base1, quality)
        }
        (false, false) => match (disagreement, quality1.cmp(&quality2)) {
            (DisagreementStrategy::MaskBoth, _) | (_, Ordering::Equal) => {
                (NO_CALL, NO_CALL_QUALITY)
            }
            (DisagreementStrategy::Consensus, Ordering::Greater) => {
                (base1, (quality1 - quality2).max(NO_CALL_QUALITY))
            }
            (DisagreementStrategy::Consensus, Ordering::Less) => {
                (base2, (quality2 - quality1).max(NO_CALL_QUALITY))
            }
            (DisagreementStrategy::MaskLowerQual, Ordering::Greater) => first,
            (DisagreementStrategy::MaskLowerQual, Ordering::Less) => second,
        },
    }
}
