//! The distance from a read's base to the template's other end, in bases of the template.
//!
//! The other end of a read's template is the 5′ end of its mate in an FR pair. The distance counts
//! the template's bases between a reference position and that end, walking the mate's CIGAR from
//! its `MC` tag and the read's own, so an indel in either read changes it by its length. Where
//! the mate holds no base, a position both reads align to serves as a landmark, and where they
//! align to none in common, the reference between them is counted one base per position.
//! Soft-clipped bases are template bases and hard-clipped bases, which the record no longer
//! holds, are not. The template length (TLEN) is never read.
//!
//! A pair is FR as htsjdk 5.0.0's `SamPairUtil.getPairOrientation` classifies it
//! (samtools/htsjdk#1771), which fgbio pins: the forward read's aligned 5′ position is at or
//! before the reverse read's.

use std::io;

use bstr::ByteSlice;
use noodles::sam;
use noodles::sam::alignment::Record;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::value::Array;
use noodles::sam::alignment::record::data::field::{Tag, Value};

use crate::error::{Error, Result};
use crate::pileup::{EntryKind, record_name};

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
    /// The alignment of CIGAR operators placed at a 0-based start, which ends before it starts
    /// when no operator consumes the reference.
    pub(crate) fn new(start: i64, ops: Vec<(Kind, i64)>) -> Self {
        let span: i64 = ops
            .iter()
            .filter(|(kind, _)| kind.consumes_reference())
            .map(|(_, len)| len)
            .sum();
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
        Self {
            start,
            end: start + span - 1,
            leading,
            trailing,
            length,
            ops,
        }
    }

    /// Whether the alignment holds a base or a deletion on the reference.
    pub(crate) fn spans_reference(&self) -> bool {
        self.end >= self.start
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

/// The operators of a CIGAR string, or `None` for none, for text that is not a CIGAR string, or
/// for an operator longer than BAM allows.
pub(crate) fn parse_cigar(text: &[u8]) -> Option<Vec<(Kind, i64)>> {
    let ops = sam::record::Cigar::new(text)
        .iter()
        .map(|op| op.ok().filter(|op| op.len() <= LONGEST_OPERATOR))
        .map(|op| op.map(|op| (op.kind(), op.len() as i64)))
        .collect::<Option<Vec<_>>>()?;
    (!ops.is_empty()).then_some(ops)
}

/// The distance of a read's base, at a query offset, from the read's 5′ end, in bases as
/// sequenced: the offset for a forward read and counted from the other end for a reverse one, so
/// soft-clipped bases count and hard-clipped ones do not. It is `None` for an offset past the read.
pub fn five_prime_distance<R: Record + ?Sized>(
    record: &R,
    query_offset: usize,
) -> Result<Option<usize>> {
    let length = record.cigar().read_length()?;
    let reverse = record.flags()?.is_reverse_complemented();
    Ok(from_five_prime(reverse, length, query_offset))
}

/// The distance from the 5′ end of a read of a query length on a strand of what an entry holds:
/// of its base, at a query offset, or for a deletion or a skip, the read's bases sequenced before
/// it, given the offset of its next base, `None` with no next base.
pub(crate) fn entry_from_five_prime(
    kind: EntryKind,
    reverse: bool,
    length: usize,
    offset: Option<usize>,
) -> Option<usize> {
    match kind {
        EntryKind::Base => from_five_prime(reverse, length, offset?),
        EntryKind::Deletion | EntryKind::Skip => {
            let next = offset.unwrap_or(length);
            if reverse {
                length.checked_sub(next)
            } else {
                (next <= length).then_some(next)
            }
        }
        EntryKind::Insertion => None,
    }
}

/// The distance of a query offset from the 5′ end of a read of a query length on a strand.
pub(crate) fn from_five_prime(reverse: bool, length: usize, offset: usize) -> Option<usize> {
    if reverse {
        length.checked_sub(offset + 1)
    } else {
        (offset < length).then_some(offset)
    }
}

/// Whether a record is a read of an FR pair, as htsjdk 5.0.0's `SamPairUtil.getPairOrientation`
/// classifies it (samtools/htsjdk#1771), which fgbio pins.
///
/// A pair is FR when its reads are paired, mapped to one contig, and on opposite strands, and the
/// forward read's 5′ position is at or before the reverse read's, so a pair whose 5′ ends coincide
/// is FR. Positions are aligned, not unclipped: a reverse record's aligned end is compared with its
/// mate's start, and a forward record's start with its mate's aligned end, from the `MC` tag and
/// never from the template length (TLEN). An end is the start plus the reference span less one,
/// so a read that spans no reference ends before it starts. The header resolves the record's
/// contigs, which a BAM record holds as indices.
///
/// A forward record of a pair otherwise FR with no `MC` tag, or one that is not a CIGAR string, is
/// an error naming the record.
pub fn is_fr_pair<R: Record + ?Sized>(record: &R, header: &sam::Header) -> Result<bool> {
    Ok(pairing(record, header)?.is_some())
}

/// A read of an FR pair: its 0-based start, its mate's, and, for a forward read, its mate's
/// alignment, read from `MC` to classify the pair.
struct Pairing {
    start: i64,
    mate_start: i64,
    mate: Option<Alignment>,
}

/// The pairing of a record of an FR pair, or `None` for any other record.
fn pairing<R: Record + ?Sized>(record: &R, header: &sam::Header) -> Result<Option<Pairing>> {
    let flags = record.flags()?;
    let reverse = flags.is_reverse_complemented();
    let (id, mate_id) = (
        record.reference_sequence_id(header).transpose()?,
        record.mate_reference_sequence_id(header).transpose()?,
    );
    if !flags.is_segmented()
        || flags.is_unmapped()
        || flags.is_mate_unmapped()
        || reverse == flags.is_mate_reverse_complemented()
        || id.is_none()
        || id != mate_id
    {
        return Ok(None);
    }
    let (Some(start), Some(mate_start)) = (
        record.alignment_start().transpose()?,
        record.mate_alignment_start().transpose()?,
    ) else {
        return Ok(None);
    };
    let (start, mate_start) = (
        usize::from(start) as i64 - 1,
        usize::from(mate_start) as i64 - 1,
    );
    let (fr, mate) = if reverse {
        let end = start + record.cigar().alignment_span()? as i64 - 1;
        (mate_start <= end, None)
    } else {
        let mate = Alignment::new(mate_start, mate_cigar(record)?);
        (start <= mate.end, Some(mate))
    };
    Ok(fr.then_some(Pairing {
        start,
        mate_start,
        mate,
    }))
}

/// The operators of a record's mate's CIGAR, from its `MC` tag.
fn mate_cigar<R: Record + ?Sized>(record: &R) -> Result<Vec<(Kind, i64)>> {
    let Some(value) = record.data().get(&Tag::MATE_CIGAR).transpose()? else {
        return Err(Error::MissingMateCigar {
            name: record_name(record),
        });
    };
    match &value {
        Value::String(text) => parse_cigar(text),
        _ => None,
    }
    .ok_or_else(|| invalid_mate_cigar(record, &value))
}

/// The error of a record whose `MC` tag holds a value that is not a usable CIGAR string.
fn invalid_mate_cigar<R: Record + ?Sized>(record: &R, value: &Value<'_>) -> Error {
    Error::InvalidMateCigar {
        name: record_name(record),
        value: value_text(value),
    }
}

/// A tag's value as SAM text spells it after the type, such as `1,-2` for an array of 1 and -2.
fn value_text(value: &Value<'_>) -> String {
    fn join<T: ToString>(values: impl Iterator<Item = io::Result<T>>) -> String {
        let texts: Vec<String> = values
            .map(|value| value.map_or_else(|_| "?".to_owned(), |value| value.to_string()))
            .collect();
        texts.join(",")
    }
    match value {
        Value::Character(character) => char::from(*character).to_string(),
        Value::Int8(number) => number.to_string(),
        Value::UInt8(number) => number.to_string(),
        Value::Int16(number) => number.to_string(),
        Value::UInt16(number) => number.to_string(),
        Value::Int32(number) => number.to_string(),
        Value::UInt32(number) => number.to_string(),
        Value::Float(number) => number.to_string(),
        Value::String(text) | Value::Hex(text) => text.to_str_lossy().into_owned(),
        Value::Array(Array::Int8(values)) => join(values.iter()),
        Value::Array(Array::UInt8(values)) => join(values.iter()),
        Value::Array(Array::Int16(values)) => join(values.iter()),
        Value::Array(Array::UInt16(values)) => join(values.iter()),
        Value::Array(Array::Int32(values)) => join(values.iter()),
        Value::Array(Array::UInt32(values)) => join(values.iter()),
        Value::Array(Array::Float(values)) => join(values.iter()),
    }
}

/// The number of the template's bases between a 0-based reference position of a record and the
/// 5′ end of its mate in an FR pair: 0 at the mate's 5′ end.
///
/// It walks the record's CIGAR and the mate's, from its `MC` tag, so an indel counts by its
/// length; soft-clipped bases count and hard-clipped bases do not, and the template length (TLEN)
/// is never read. It is `None` for any record [`is_fr_pair`] does not call a read of an FR pair,
/// such as a fragment, a read whose mate is unmapped or on another contig, or a read of an
/// outward-facing pair, and for a read of an FR pair only at a position past the mate's 5′ end.
/// The header resolves the record's contigs, which a BAM record holds as indices.
///
/// A record of an FR pair with no `MC` tag, or one that is not a CIGAR string spanning at least
/// one base, is an error naming the record.
pub fn template_end_distance<R: Record + ?Sized>(
    record: &R,
    header: &sam::Header,
    position: usize,
) -> Result<Option<usize>> {
    Ok(ends(record, header)?.and_then(|ends| ends.distance(position)))
}

/// A read of an FR pair and its mate laid along the reference, from which the template's bases
/// between a position and the mate's 5′ end are counted.
#[derive(Clone, Debug)]
pub(crate) struct Ends {
    read: Alignment,
    mate: Alignment,
    reverse: bool,
}

impl Ends {
    /// The template's bases between a 0-based position and the mate's 5′ end, or `None` past it.
    pub(crate) fn distance(&self, position: usize) -> Option<usize> {
        usize::try_from(distance(
            &self.read,
            &self.mate,
            self.reverse,
            position as i64,
        ))
        .ok()
    }
}

/// The ends of a read of an FR pair, as [`template_end_distance`] counts from them, or `None` for
/// any other record or a read that spans no reference.
pub(crate) fn ends<R: Record + ?Sized>(record: &R, header: &sam::Header) -> Result<Option<Ends>> {
    let Some(pairing) = pairing(record, header)? else {
        return Ok(None);
    };
    let mate = match pairing.mate {
        Some(mate) => mate,
        None => Alignment::new(pairing.mate_start, mate_cigar(record)?),
    };
    if !mate.spans_reference() {
        let value = record
            .data()
            .get(&Tag::MATE_CIGAR)
            .transpose()?
            .unwrap_or(Value::String(b"".as_bstr()));
        return Err(invalid_mate_cigar(record, &value));
    }
    let ops = record
        .cigar()
        .iter()
        .map(|op| op.map(|op| (op.kind(), op.len() as i64)))
        .collect::<io::Result<Vec<_>>>()?;
    let read = Alignment::new(pairing.start, ops);
    let reverse = record.flags()?.is_reverse_complemented();
    Ok(read.spans_reference().then_some(Ends {
        read,
        mate,
        reverse,
    }))
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
