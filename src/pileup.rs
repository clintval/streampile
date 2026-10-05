use bstr::BStr;
use noodles::bam;
use noodles::sam::alignment::record::Flags;

use crate::auxiliary::{self, AuxElement, AuxValue, Field};
use crate::error::Result;
use crate::footprint::Footprint;

/// The base quality of every base of a read with no stored qualities (QUAL `*`), as in htslib.
pub const MISSING_BASE_QUALITY: u8 = 255;

pub(crate) const NONE: u32 = u32::MAX;

/// What a read holds at a pileup's position.
///
/// A [`Skip`](EntryKind::Skip) is a reference skip, the CIGAR `N` operator, not an `N` base,
/// which is a [`Base`](EntryKind::Base).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum EntryKind {
    /// An aligned base, a match or a mismatch.
    Base,
    /// A deletion.
    Deletion,
    /// Bases inserted right after the position.
    Insertion,
    /// A reference skip.
    Skip,
}

impl EntryKind {
    /// The kind's name as the Python package spells it.
    pub fn as_str(self) -> &'static str {
        match self {
            EntryKind::Base => "base",
            EntryKind::Deletion => "deletion",
            EntryKind::Insertion => "insertion",
            EntryKind::Skip => "skip",
        }
    }
}

/// One record held by a builder, with what was worked out once when it was read.
#[derive(Debug)]
pub(crate) struct LiveRecord {
    pub record: bam::Record,
    pub reference_id: usize,
    pub start: i64,
    pub end: i64,
    pub flags: Flags,
    pub footprint: Footprint,
    pub other_end: Option<i64>,
    pub fields: Vec<Option<Field>>,
    pub template: u32,
}

impl Default for LiveRecord {
    fn default() -> Self {
        Self {
            record: bam::Record::default(),
            reference_id: usize::MAX,
            start: -1,
            end: 0,
            flags: Flags::empty(),
            footprint: Footprint::default(),
            other_end: None,
            fields: Vec::new(),
            template: NONE,
        }
    }
}

/// An entry as a builder stores it: the record's slot and what it holds at the position.
#[derive(Clone, Copy, Debug)]
pub(crate) struct RawEntry {
    pub slot: u32,
    pub kind: EntryKind,
    pub offset: u32,
    pub length: u32,
}

/// The reads at one reference position, borrowed from the builder that built it.
///
/// Entries come in input order. A read with an insertion right after the position appears again
/// as an insertion entry, after its entry at the position, so one read can have two entries. A
/// read that opens with an insertion has one at the position before its first aligned base.
///
/// Every entry of an accepted read is kept: `min_base_quality` applies only to
/// [`filtered_depth`](Pileup::filtered_depth), [`bases`](Pileup::bases), and
/// [`qualities`](Pileup::qualities), and to choosing the mate a builder keeps when it leaves out
/// overlapping mates.
#[derive(Clone, Copy, Debug)]
pub struct Pileup<'a> {
    pub(crate) reference_sequence_id: usize,
    pub(crate) reference_sequence_name: &'a BStr,
    pub(crate) position: usize,
    pub(crate) min_base_quality: u8,
    pub(crate) entries: &'a [RawEntry],
    pub(crate) slots: &'a [LiveRecord],
    pub(crate) aux_tags: &'a [[u8; 2]],
}

impl<'a> Pileup<'a> {
    /// The index of the contig in the header.
    pub fn reference_sequence_id(&self) -> usize {
        self.reference_sequence_id
    }

    /// The name of the contig.
    pub fn reference_sequence_name(&self) -> &'a BStr {
        self.reference_sequence_name
    }

    /// The 0-based position on the contig.
    pub fn position(&self) -> usize {
        self.position
    }

    /// The base quality below which bases are left out of the filtered views.
    pub fn min_base_quality(&self) -> u8 {
        self.min_base_quality
    }

    /// The number of entries, insertion entries included.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Whether the pileup has no entries.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// The entry at an index.
    pub fn get(&self, index: usize) -> Option<PileupEntry<'a>> {
        self.entries.get(index).map(|raw| self.entry(*raw))
    }

    /// Every entry, in input order.
    pub fn iter(&self) -> impl ExactSizeIterator<Item = PileupEntry<'a>> + 'a {
        let pileup = *self;
        self.entries.iter().map(move |raw| pileup.entry(*raw))
    }

    fn entry(&self, raw: RawEntry) -> PileupEntry<'a> {
        PileupEntry {
            live: &self.slots[raw.slot as usize],
            raw,
            position: self.position as i64,
            aux_tags: self.aux_tags,
        }
    }

    /// The number of reads with a base, a deletion, or a skip here.
    ///
    /// Quality is ignored, and insertion entries are not counted, as in htslib's column depth.
    pub fn unfiltered_depth(&self) -> usize {
        self.entries
            .iter()
            .filter(|raw| raw.kind != EntryKind::Insertion)
            .count()
    }

    /// The number of reads with a base or a deletion here at the quality floor.
    ///
    /// A deletion is judged by the quality of the read's next base. A skip has no quality, so it
    /// is never counted.
    pub fn filtered_depth(&self) -> usize {
        let floor = self.min_base_quality;
        self.iter()
            .filter(|entry| entry.kind() != EntryKind::Insertion && entry.passes(floor))
            .count()
    }

    /// The upper-cased bases here at the quality floor, in entry order.
    pub fn bases(&self) -> impl Iterator<Item = u8> + 'a {
        let floor = self.min_base_quality;
        self.iter()
            .filter(move |entry| entry.passes(floor))
            .filter_map(|entry| entry.base())
    }

    /// The qualities of the bases here at the quality floor, lined up with
    /// [`bases`](Pileup::bases).
    pub fn qualities(&self) -> impl Iterator<Item = u8> + 'a {
        let floor = self.min_base_quality;
        self.iter()
            .filter(|entry| entry.kind() == EntryKind::Base)
            .filter_map(|entry| entry.quality())
            .filter(move |&quality| quality >= floor)
    }
}

/// One read at one pileup position.
#[derive(Clone, Copy, Debug)]
pub struct PileupEntry<'a> {
    live: &'a LiveRecord,
    raw: RawEntry,
    position: i64,
    aux_tags: &'a [[u8; 2]],
}

impl<'a> PileupEntry<'a> {
    /// The record.
    pub fn record(&self) -> &'a bam::Record {
        &self.live.record
    }

    /// The record's flags.
    pub fn flags(&self) -> Flags {
        self.live.flags
    }

    /// What the read holds here.
    pub fn kind(&self) -> EntryKind {
        self.raw.kind
    }

    /// Whether the read has a deletion here.
    pub fn is_deletion(&self) -> bool {
        self.raw.kind == EntryKind::Deletion
    }

    /// Whether this is an insertion right after the position.
    pub fn is_insertion(&self) -> bool {
        self.raw.kind == EntryKind::Insertion
    }

    /// Whether the read skips over the position with the CIGAR `N` operator.
    pub fn is_skip(&self) -> bool {
        self.raw.kind == EntryKind::Skip
    }

    /// Whether the read holds an `N` base, a no-call, here.
    pub fn is_no_call(&self) -> bool {
        self.base() == Some(b'N')
    }

    /// Whether the read is aligned to the reverse strand.
    pub fn is_reverse(&self) -> bool {
        self.live.flags.is_reverse_complemented()
    }

    /// The 0-based query offset of the read's base here, for a base entry.
    pub fn query_position(&self) -> Option<usize> {
        (self.raw.kind == EntryKind::Base).then_some(self.raw.offset as usize)
    }

    /// The query offset of the read's base here, or for a deletion of the read's next base, as
    /// htslib reports it, or `None` for a skip, an insertion, or a deletion no base follows.
    pub fn query_position_or_next(&self) -> Option<usize> {
        match self.raw.kind {
            EntryKind::Base | EntryKind::Deletion if self.raw.offset != NONE => {
                Some(self.raw.offset as usize)
            }
            _ => None,
        }
    }

    /// The query offset of the first inserted base of an insertion entry.
    pub fn insertion_offset(&self) -> Option<usize> {
        (self.raw.kind == EntryKind::Insertion).then_some(self.raw.offset as usize)
    }

    /// The number of inserted bases of an insertion entry, or 0 for any other.
    pub fn insertion_length(&self) -> usize {
        if self.raw.kind == EntryKind::Insertion {
            self.raw.length as usize
        } else {
            0
        }
    }

    /// The upper-cased base here, or `None` without one, as for a read with no stored bases.
    pub fn base(&self) -> Option<u8> {
        self.query_position()
            .and_then(|offset| self.live.record.sequence().get(offset))
    }

    /// The base as the sequencer read it: the complement of [`base`](PileupEntry::base) for a
    /// reverse read.
    pub fn sequenced_base(&self) -> Option<u8> {
        let base = self.base()?;
        Some(if self.is_reverse() {
            complement(base)
        } else {
            base
        })
    }

    /// The base quality here, or of the read's next base for a deletion.
    ///
    /// A read with bases but no stored qualities (QUAL `*`) has quality 255 at every base, as in
    /// htslib. It is `None` for an insertion, a skip, a deletion no base follows, or a read with
    /// no stored bases (SEQ `*`).
    pub fn quality(&self) -> Option<u8> {
        quality_of(&self.live.record, self.raw.kind, self.raw.offset)
    }

    /// Whether the entry has a quality at the floor: a base or a deletion followed by a base.
    pub fn passes(&self, min_base_quality: u8) -> bool {
        self.quality()
            .is_some_and(|quality| quality >= min_base_quality)
    }

    /// The upper-cased inserted bases of an insertion entry, or `None` for any other entry or a
    /// read with no stored bases.
    pub fn inserted_bases(&self) -> Option<impl Iterator<Item = u8> + 'a> {
        let offset = self.insertion_offset()?;
        let sequence = self.live.record.sequence();
        if sequence.is_empty() {
            return None;
        }
        let end = (offset + self.insertion_length()).min(sequence.len());
        Some((offset.min(end)..end).filter_map(move |index| sequence.get(index)))
    }

    /// The qualities of the inserted bases of an insertion entry: 255 for a read with no stored
    /// qualities, and `None` for any other entry or a read with no stored bases.
    pub fn inserted_qualities(&self) -> Option<impl Iterator<Item = u8> + 'a> {
        let offset = self.insertion_offset()?;
        let record = &self.live.record;
        let stored = record.sequence().len();
        if stored == 0 {
            return None;
        }
        let qualities = record.quality_scores().as_bytes();
        let length = self.insertion_length();
        let (bytes, missing): (&'a [u8], usize) = if qualities.is_empty() {
            (&[], length)
        } else {
            let end = (offset + length).min(qualities.len());
            (&qualities[offset.min(end)..end], 0)
        };
        Some(
            bytes
                .iter()
                .copied()
                .chain(std::iter::repeat_n(MISSING_BASE_QUALITY, missing)),
        )
    }

    /// The distance of this base from the read's 5′ end, in bases as sequenced.
    ///
    /// It is the query offset for a forward read, counted from the other end for a reverse read,
    /// so soft-clipped bases count, and 0 is the first base sequenced: fgbio's
    /// `positionInReadInReadOrder` minus one. It is `None` for an entry with no base.
    pub fn five_prime_distance(&self) -> Option<usize> {
        let offset = self.query_position()?;
        if self.is_reverse() {
            (self.live.footprint.query_length as usize).checked_sub(offset + 1)
        } else {
            Some(offset)
        }
    }

    /// The distance on the reference from this position to the template's other end, the 5′ end
    /// of the mate of a read in an FR pair: 0 at the mate's 5′ end.
    ///
    /// For a reverse read, the mate's 5′ end is its alignment start. For a forward read, it is the
    /// end of the mate's alignment, from its start and its `MC` tag, or with no `MC` tag, from
    /// the read's own start and its template length (TLEN). It is `None` for a fragment, a read
    /// whose mate is unmapped or on another contig, a pair that is not FR, and a position past
    /// the mate's 5′ end, where a read runs through its mate.
    pub fn template_end_distance(&self) -> Option<usize> {
        let other = self.live.other_end?;
        let position = self.position;
        let distance = if self.is_reverse() {
            position - other
        } else {
            other - position
        };
        usize::try_from(distance).ok()
    }

    /// The value of one of the record's auxiliary fields, borrowed from the record.
    ///
    /// A tag the builder was asked to index with
    /// [`index_aux_tags`](crate::StreamingPileupBuilder::index_aux_tags) is found without a
    /// search; any other is searched for in the record's fields.
    pub fn aux(&self, tag: [u8; 2]) -> Result<Option<AuxValue<'a>>> {
        let data = self.live.record.data().as_bytes();
        match self.aux_tags.iter().position(|indexed| *indexed == tag) {
            Some(index) => match self.live.fields.get(index).copied().flatten() {
                Some(field) => field.value(data).map(Some),
                None => Ok(None),
            },
            None => auxiliary::find(data, tag),
        }
    }

    /// The element at this base's query offset of a per-base string or array field, such as a
    /// per-base depth, or `None` for an entry with no base or a record without the field.
    pub fn aux_at(&self, tag: [u8; 2]) -> Result<Option<AuxElement>> {
        let Some(offset) = self.query_position() else {
            return Ok(None);
        };
        Ok(self.aux(tag)?.and_then(|value| value.get(offset)))
    }
}

/// The quality of the base at a query offset of a base or deletion entry, as
/// [`PileupEntry::quality`] reports it.
pub(crate) fn quality_of(record: &bam::Record, kind: EntryKind, offset: u32) -> Option<u8> {
    if offset == NONE || !matches!(kind, EntryKind::Base | EntryKind::Deletion) {
        return None;
    }
    let offset = offset as usize;
    let qualities = record.quality_scores().as_bytes();
    if qualities.is_empty() {
        (offset < record.sequence().len()).then_some(MISSING_BASE_QUALITY)
    } else {
        qualities.get(offset).copied()
    }
}

fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'M' => b'K',
        b'K' => b'M',
        b'R' => b'Y',
        b'Y' => b'R',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        other => other,
    }
}
