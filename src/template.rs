//! The distance from a read's base to the template's other end, in bases of the template.
//!
//! The other end of a read's template is the 5′ end of its mate in an FR pair. The distance counts
//! the template's bases between a reference position and that end, walking the mate's CIGAR from
//! its `MC` tag and the read's own, so an indel in either read changes it by its length. Where
//! the mate holds no base, a position both reads align to serves as a landmark, and where they
//! align to none in common, the reference between them is counted one base per position.
//! Soft-clipped bases are template bases and hard-clipped bases, which the record no longer
//! holds, are not. The template length (TLEN) is never read.

use noodles::bam;
use noodles::sam;
use noodles::sam::alignment::record::cigar::op::Kind;

use crate::auxiliary::{self, AuxValue};
use crate::error::{Error, Result};
use crate::pileup::record_name;

const LONGEST_OPERATOR: usize = (1 << 28) - 1;

/// One alignment's query bases laid along the reference.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Alignment {
    start: i64,
    end: i64,
    leading: i64,
    trailing: i64,
    length: i64,
    ops: Vec<(Kind, i64)>,
}

impl Alignment {
    /// The alignment of CIGAR operators placed at a 0-based start, or `None` without a base on
    /// the reference.
    pub(crate) fn new(start: i64, ops: Vec<(Kind, i64)>) -> Option<Self> {
        let span: i64 = ops
            .iter()
            .filter(|(kind, _)| kind.consumes_reference())
            .map(|(_, len)| len)
            .sum();
        if span == 0 {
            return None;
        }
        let soft = |op: &&(Kind, i64)| op.0 == Kind::SoftClip;
        let clip = |op: &&(Kind, i64)| matches!(op.0, Kind::SoftClip | Kind::HardClip);
        let leading = ops
            .iter()
            .take_while(clip)
            .filter(soft)
            .map(|op| op.1)
            .sum();
        let trailing = ops
            .iter()
            .rev()
            .take_while(clip)
            .filter(soft)
            .map(|op| op.1)
            .sum();
        let length = ops
            .iter()
            .filter(|(kind, _)| kind.consumes_read())
            .map(|(_, len)| len)
            .sum();
        Some(Self {
            start,
            end: start + span - 1,
            leading,
            trailing,
            length,
            ops,
        })
    }

    /// The alignment an `MC` value describes for a mate at a 0-based start, or `None` for a
    /// value that is not a CIGAR string spanning at least one base with operators no longer than
    /// BAM allows.
    pub(crate) fn of_mate(start: i64, value: AuxValue<'_>) -> Option<Self> {
        let AuxValue::String(text) = value else {
            return None;
        };
        let ops = sam::record::Cigar::new(text)
            .iter()
            .map(|op| op.ok().filter(|op| op.len() <= LONGEST_OPERATOR))
            .map(|op| op.map(|op| (op.kind(), op.len() as i64)))
            .collect::<Option<Vec<_>>>()?;
        Self::new(start, ops)
    }

    /// The query offset at a 0-based reference position, and whether a base lies there.
    ///
    /// Inside a deletion or a skip it is the offset of the next base, with no base. Before the
    /// alignment and after it, leading and trailing soft clips, and then the reference, are
    /// counted one base per position, so offsets run below 0 and past the query.
    pub(crate) fn offset(&self, position: i64) -> (i64, bool) {
        if position < self.start {
            return (self.leading - (self.start - position), true);
        }
        if position > self.end {
            return (
                self.length - self.trailing + (position - self.end - 1),
                true,
            );
        }
        let (mut query, mut reference) = (self.leading, self.start);
        for &(kind, len) in self
            .ops
            .iter()
            .skip_while(|(kind, _)| matches!(kind, Kind::SoftClip | Kind::HardClip))
        {
            match kind {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    if position < reference + len {
                        return (query + position - reference, true);
                    }
                    query += len;
                    reference += len;
                }
                Kind::Insertion => query += len,
                Kind::Deletion | Kind::Skip => {
                    if position < reference + len {
                        return (query, false);
                    }
                    reference += len;
                }
                Kind::SoftClip | Kind::HardClip | Kind::Pad => {}
            }
        }
        (query, false)
    }

    /// The query bases strictly before a reference position.
    fn before(&self, position: i64) -> i64 {
        self.offset(position).0
    }

    /// The query bases strictly after a reference position.
    fn after(&self, position: i64) -> i64 {
        let (offset, base) = self.offset(position);
        self.length - offset - i64::from(base)
    }

    /// The query bases strictly between two reference positions, `from` before `to`.
    fn between(&self, from: i64, to: i64) -> i64 {
        let (from, base) = self.offset(from);
        self.offset(to).0 - from - i64::from(base)
    }
}

/// The template's bases between a 0-based reference position of a record and the 5′ end of its
/// mate in an FR pair, or `None` for a fragment, a read whose mate is unmapped or on another
/// contig, a pair that is not FR, and a position past the mate's 5′ end.
///
/// A record of an FR pair without an `MC` tag, or with one that is not a usable CIGAR string, is
/// an error naming the record.
pub(crate) fn template_end_distance(record: &bam::Record, position: i64) -> Result<Option<usize>> {
    let flags = record.flags();
    let reverse = flags.is_reverse_complemented();
    if !flags.is_segmented()
        || flags.is_unmapped()
        || flags.is_mate_unmapped()
        || reverse == flags.is_mate_reverse_complemented()
        || record.reference_sequence_id().transpose()?
            != record.mate_reference_sequence_id().transpose()?
    {
        return Ok(None);
    }
    let (Some(start), Some(mate_start)) = (
        record.alignment_start().transpose()?,
        record.mate_alignment_start().transpose()?,
    ) else {
        return Ok(None);
    };
    let ops = record
        .cigar()
        .iter()
        .map(|op| op.map(|op| (op.kind(), op.len() as i64)))
        .collect::<std::io::Result<Vec<_>>>()?;
    let Some(read) = Alignment::new(usize::from(start) as i64 - 1, ops) else {
        return Ok(None);
    };
    let name = || record_name(record);
    let Some(value) = auxiliary::find(record.data().as_bytes(), *b"MC")? else {
        return Err(Error::MissingMateCigar { name: name() });
    };
    let Some(mate) = Alignment::of_mate(usize::from(mate_start) as i64 - 1, value) else {
        return Err(Error::InvalidMateCigar {
            name: name(),
            value: match value {
                AuxValue::String(text) => String::from_utf8_lossy(text).into_owned(),
                other => format!("{other:?}"),
            },
        });
    };
    Ok(usize::try_from(distance(&read, &mate, reverse, position)).ok())
}

/// The template's bases between a position of a read and its mate's 5′ end, which is past the
/// mate's end for a forward read and before its start for a reverse one; negative past that end.
pub(crate) fn distance(read: &Alignment, mate: &Alignment, reverse: bool, position: i64) -> i64 {
    let (_, base) = read.offset(position);
    let base = i64::from(base);
    if !reverse {
        if position < mate.start {
            read.between(position, mate.start) + mate.length - mate.leading
        } else if position <= mate.end {
            mate.after(position)
        } else {
            mate.after(mate.end) - read.between(mate.end, position) - base
        }
    } else if position > mate.end {
        mate.before(mate.end) + 1 + read.between(mate.end, position)
    } else if position >= mate.start {
        mate.before(position)
    } else {
        mate.leading - read.between(position, mate.start) - base
    }
}
