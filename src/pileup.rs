use std::collections::HashMap;
use std::hash::{BuildHasherDefault, DefaultHasher, Hash, Hasher};
use std::sync::OnceLock;

use bstr::{BStr, ByteSlice};
use noodles::bam;
use noodles::sam::alignment::record::Flags;

use crate::auxiliary::{self, AuxElement, AuxValue, Field};
use crate::error::Result;
use crate::footprint::Footprint;
use crate::overlap::{AgreementStrategy, DisagreementStrategy, Observation, Vote};
use crate::source::AlignmentRecord;
use crate::template::{self, Ends};

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
pub(crate) struct LiveRecord<R = bam::Record> {
    pub record: R,
    pub reference_id: usize,
    pub start: i64,
    pub end: i64,
    pub flags: Flags,
    pub footprint: Footprint,
    pub fields: Vec<Option<Field>>,
    pub derived: Derived,
}

impl<R: Default> Default for LiveRecord<R> {
    fn default() -> Self {
        Self {
            record: R::default(),
            reference_id: usize::MAX,
            start: -1,
            end: 0,
            flags: Flags::empty(),
            footprint: Footprint::default(),
            fields: Vec::new(),
            derived: Derived::default(),
        }
    }
}

impl<R: AlignmentRecord> LiveRecord<R> {
    /// The template of the record held in a slot.
    pub fn template(&self, slot: u32) -> TemplateName<'_> {
        TemplateName::of(self.record.bam(), &self.derived, slot as usize)
    }
}

/// What is worked out about a record only when it is first asked for, and then kept.
#[derive(Debug, Default)]
pub(crate) struct Derived {
    name_hash: OnceLock<u64>,
    ends: OnceLock<Option<Ends>>,
}

impl Derived {
    /// The hash of the record's name.
    fn name_hash(&self, record: &bam::Record) -> u64 {
        *self.name_hash.get_or_init(|| {
            let mut hasher = DefaultHasher::new();
            hasher.write(name_of(record));
            hasher.finish()
        })
    }

    /// The record's ends, kept once worked out, or the error of a record missing them, worked
    /// out again each time.
    pub(crate) fn ends(&self, record: &bam::Record) -> Result<Option<&Ends>> {
        if let Some(ends) = self.ends.get() {
            return Ok(ends.as_ref());
        }
        let ends = template::ends(record, &noodles::sam::Header::default())?;
        Ok(self.ends.get_or_init(|| ends).as_ref())
    }

    /// Whether the record is a read of an FR pair, from its ends where it has them.
    pub(crate) fn is_fr_pair(&self, record: &bam::Record) -> Result<bool> {
        match self.ends(record) {
            Ok(ends) => Ok(ends.is_some()),
            Err(_) => crate::is_fr_pair(record, &noodles::sam::Header::default()),
        }
    }
}

/// A read's template: its name and a hash of it worked out once, or for a read with no name,
/// which is a template of its own, a number unique to the read.
#[derive(Clone, Copy, Debug)]
pub(crate) struct TemplateName<'n> {
    hash: u64,
    name: Option<&'n [u8]>,
    read: usize,
}

impl<'n> TemplateName<'n> {
    /// The template of a record, and of a read with no name, a number unique to the read.
    pub(crate) fn of(record: &'n bam::Record, derived: &Derived, read: usize) -> Self {
        let name = record.name().map(|name| name.as_bytes());
        Self {
            hash: match name {
                Some(_) => derived.name_hash(record),
                None => read as u64,
            },
            name,
            read,
        }
    }
}

impl Hash for TemplateName<'_> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        state.write_u64(self.hash);
    }
}

impl PartialEq for TemplateName<'_> {
    fn eq(&self, other: &Self) -> bool {
        match (self.name, other.name) {
            (Some(name), Some(other)) => name == other,
            (None, None) => self.read == other.read,
            _ => false,
        }
    }
}

impl Eq for TemplateName<'_> {}

/// A hasher of [`TemplateName`] keys, which are hashed already.
#[derive(Default)]
struct Prehashed(u64);

impl Hasher for Prehashed {
    fn finish(&self) -> u64 {
        self.0
    }

    fn write(&mut self, bytes: &[u8]) {
        for &byte in bytes {
            self.0 = self.0.rotate_left(8) ^ u64::from(byte);
        }
    }

    fn write_u64(&mut self, hash: u64) {
        self.0 = hash;
    }
}

/// A record's name, `*` for a record with none.
pub(crate) fn name_of(record: &bam::Record) -> &[u8] {
    record.name().map_or(b"*", |name| name.as_bytes())
}

/// The template of each entry, numbered in the order of each template's first entry, or `None`
/// for an entry that `name` gives no template.
pub(crate) fn number_templates<'n, T>(
    entries: &'n [T],
    name: impl Fn(&'n T) -> Option<TemplateName<'n>>,
) -> Vec<Option<usize>> {
    let mut numbers: HashMap<TemplateName<'n>, usize, BuildHasherDefault<Prehashed>> =
        HashMap::with_capacity_and_hasher(entries.len(), BuildHasherDefault::default());
    entries
        .iter()
        .map(|entry| {
            let name = name(entry)?;
            let next = numbers.len();
            Some(*numbers.entry(name).or_insert(next))
        })
        .collect()
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
/// [`qualities`](Pileup::qualities). [`templates`](Pileup::templates) sees each template once,
/// calling the bases of overlapping mates into one.
///
/// `R` is the type of record the builder's source reads: [`bam::Record`] for a BAM reader.
#[derive(Debug)]
pub struct Pileup<'a, R = bam::Record> {
    pub(crate) reference_sequence_id: usize,
    pub(crate) reference_sequence_name: &'a BStr,
    pub(crate) position: usize,
    pub(crate) min_base_quality: u8,
    pub(crate) entries: &'a [RawEntry],
    pub(crate) slots: &'a [LiveRecord<R>],
    pub(crate) aux_tags: &'a [[u8; 2]],
}

impl<R> Clone for Pileup<'_, R> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<R> Copy for Pileup<'_, R> {}

impl<'a, R: AlignmentRecord> Pileup<'a, R> {
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
    pub fn get(&self, index: usize) -> Option<PileupEntry<'a, R>> {
        self.entries.get(index).map(|raw| self.entry(*raw))
    }

    /// Every entry, in input order.
    pub fn iter(&self) -> impl ExactSizeIterator<Item = PileupEntry<'a, R>> + 'a {
        let pileup = *self;
        self.entries.iter().map(move |raw| pileup.entry(*raw))
    }

    fn entry(&self, raw: RawEntry) -> PileupEntry<'a, R> {
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

    /// One observation per template here, its reads grouped by name, in the order of each
    /// template's first entry. A read with no name is a template of its own.
    ///
    /// Where two reads of a template hold bases, they are called into one as fgbio's
    /// `CallOverlappingConsensusBases` calls them, as fgumi implements it: `agreement` makes the
    /// quality of equal bases and `disagreement` the base and quality of different ones. Where the
    /// strategies leave a read's base or quality unchanged, the template takes the higher quality.
    ///
    /// fgumi defines no more than that, so at each position:
    ///
    /// - A no-call (`N`) is left alone, as fgumi leaves it, so the other read's base stands at its
    ///   own quality, and two no-calls are an `N` at the higher quality. A no-call is still a
    ///   base, so it stands over the other read's deletion, which fgumi leaves alone too.
    /// - A read with a deletion or a skip holds no base, so a template whose other read holds one
    ///   has that base at its own quality. A template with no base is a deletion if either read
    ///   holds one, at the higher of their qualities, or else a skip.
    /// - Insertion entries are no part of a template, so a read whose only entry here is an
    ///   insertion adds nothing to its template, and is in none without another read here.
    /// - A read under the pileup's quality floor does not vote, as [`bases`](Pileup::bases),
    ///   [`qualities`](Pileup::qualities), and [`filtered_depth`](Pileup::filtered_depth) leave it
    ///   out: a base, or a deletion judged by its next base, under `min_base_quality`. So a mate
    ///   under the floor neither masks nor lowers the other mate's base. A template none of whose
    ///   reads votes holds no base or quality, and [`passes`](PileupTemplate::passes) compares a
    ///   template's quality to another floor.
    /// - A template with more than two reads here, as when supplementary records are piled up,
    ///   calls them in input order.
    pub fn templates(
        &self,
        agreement: AgreementStrategy,
        disagreement: DisagreementStrategy,
    ) -> Vec<PileupTemplate<'a, R>> {
        let slots = self.slots;
        let numbers = number_templates(self.entries, |raw| {
            (raw.kind != EntryKind::Insertion).then(|| slots[raw.slot as usize].template(raw.slot))
        });
        let entries = self
            .entries
            .iter()
            .zip(numbers)
            .filter_map(|(raw, number)| {
                let entry = self.entry(*raw);
                number.map(|number| (number, entry, entry.observation()))
            });
        templates(
            entries,
            agreement,
            disagreement,
            i64::from(self.min_base_quality),
        )
        .into_iter()
        .map(PileupTemplate)
        .collect()
    }
}

/// One read at one pileup position.
#[derive(Debug)]
pub struct PileupEntry<'a, R = bam::Record> {
    pub(crate) live: &'a LiveRecord<R>,
    raw: RawEntry,
    position: i64,
    aux_tags: &'a [[u8; 2]],
}

impl<R> Clone for PileupEntry<'_, R> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<R> Copy for PileupEntry<'_, R> {}

impl<'a, R: AlignmentRecord> PileupEntry<'a, R> {
    /// The BAM record.
    pub fn record(&self) -> &'a bam::Record {
        self.live.record.bam()
    }

    /// The record as the builder's source read it, which for a BAM reader is the BAM record.
    pub fn source_record(&self) -> &'a R {
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
            .and_then(|offset| self.record().sequence().get(offset))
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
        quality_of(self.record(), self.raw.kind, self.raw.offset)
    }

    /// Whether the entry has a quality at the floor: a base or a deletion followed by a base.
    pub fn passes(&self, min_base_quality: u8) -> bool {
        self.quality()
            .is_some_and(|quality| quality >= min_base_quality)
    }

    /// What the read holds here, as its template's vote counts it.
    pub(crate) fn observation(&self) -> Observation {
        Observation {
            kind: self.kind(),
            base: self.base(),
            quality: self.quality(),
        }
    }

    /// The upper-cased inserted bases of an insertion entry, or `None` for any other entry or a
    /// read with no stored bases.
    pub fn inserted_bases(&self) -> Option<impl Iterator<Item = u8> + 'a> {
        let offset = self.insertion_offset()?;
        let sequence = self.record().sequence();
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
        let record = self.record();
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
    /// `positionInReadInReadOrder` minus one. For a deletion or a skip, which holds no base, it is
    /// the number of the read's bases sequenced before the position, as
    /// [`template_end_distance`](PileupEntry::template_end_distance) counts the bases after it.
    /// It is `None` for an insertion entry.
    pub fn five_prime_distance(&self) -> Option<usize> {
        let length = self.live.footprint.query_length as usize;
        let offset = (self.raw.offset != NONE).then_some(self.raw.offset as usize);
        crate::template::entry_from_five_prime(self.raw.kind, self.is_reverse(), length, offset)
    }

    /// The number of the template's bases between this position and the template's other end,
    /// the 5′ end of the mate of a read in an FR pair: 0 at the mate's 5′ end.
    ///
    /// It walks the read's CIGAR and the mate's, from its `MC` tag, so an indel counts by its
    /// length; soft-clipped bases count and hard-clipped bases, absent from the records, do not.
    /// Where the mate has no base, a position both reads align to carries the count, or else the
    /// reference between them; the template length (TLEN) is never read. It is `None` for a read
    /// that [`is_fr_pair`](PileupEntry::is_fr_pair) does not call a read of an FR pair, and for a
    /// read of an FR pair only at a position past the mate's 5′ end, where a read runs through its
    /// mate.
    ///
    /// A read of an FR pair with no `MC` tag, or one that is not a CIGAR string spanning at least
    /// one base, is an error naming the read.
    pub fn template_end_distance(&self) -> Result<Option<usize>> {
        let ends = self.live.derived.ends(self.record())?;
        Ok(ends.and_then(|ends| ends.distance(self.position as usize)))
    }

    /// Whether the read is a read of an FR pair, as htsjdk 5.0.0's
    /// `SamPairUtil.getPairOrientation` classifies it, with the forward read's aligned 5′
    /// position at or before the reverse read's, as [`is_fr_pair`](crate::is_fr_pair) says.
    ///
    /// A forward read of a pair otherwise FR with no `MC` tag, or one that is not a CIGAR string,
    /// is an error naming the read.
    pub fn is_fr_pair(&self) -> Result<bool> {
        self.live.derived.is_fr_pair(self.record())
    }

    /// The value of one of the record's auxiliary fields, borrowed from the record.
    ///
    /// A tag the builder was asked to index with
    /// [`index_aux_tags`](crate::StreamingPileupBuilder::index_aux_tags) is found without a
    /// search; any other is searched for in the record's fields.
    pub fn aux(&self, tag: [u8; 2]) -> Result<Option<AuxValue<'a>>> {
        let data = self.record().data().as_bytes();
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

/// What a template needs of each of its reads.
pub(crate) trait TemplateRead {
    fn flags(&self) -> Flags;

    fn five_prime_distance(&self) -> Option<usize>;

    fn template_end_distance(&self) -> Result<Option<usize>>;

    fn is_fr_pair(&self) -> Result<bool>;
}

impl<R: AlignmentRecord> TemplateRead for PileupEntry<'_, R> {
    fn flags(&self) -> Flags {
        PileupEntry::flags(self)
    }

    fn five_prime_distance(&self) -> Option<usize> {
        PileupEntry::five_prime_distance(self)
    }

    fn template_end_distance(&self) -> Result<Option<usize>> {
        PileupEntry::template_end_distance(self)
    }

    fn is_fr_pair(&self) -> Result<bool> {
        PileupEntry::is_fr_pair(self)
    }
}

/// The reads of one template at one position, in input order, and the vote of their bases.
#[derive(Clone, Debug)]
pub(crate) struct Template<E> {
    first: E,
    second: Option<E>,
    others: Vec<E>,
    vote: Vote,
}

impl<E> Template<E> {
    fn of(entry: E) -> Self {
        Self {
            first: entry,
            second: None,
            others: Vec::new(),
            vote: Vote::default(),
        }
    }

    fn push(&mut self, entry: E) {
        if self.second.is_none() {
            self.second = Some(entry);
        } else {
            self.others.push(entry);
        }
    }

    /// The template's first entry.
    pub(crate) fn first(&self) -> &E {
        &self.first
    }

    /// The template's entries, in input order.
    pub(crate) fn entries(&self) -> impl Iterator<Item = &E> {
        std::iter::once(&self.first)
            .chain(&self.second)
            .chain(&self.others)
    }

    /// What the template holds: its kind, and the base and quality of its reads' vote.
    pub(crate) fn called(&self) -> Observation {
        self.vote.result()
    }
}

impl<E: TemplateRead> Template<E> {
    /// Whether the template's first read is aligned to the reverse strand, read from the second
    /// read's flags without the first.
    pub(crate) fn is_reverse(&self) -> bool {
        match self.first_read() {
            Some(first) => first.flags().is_reverse_complemented(),
            None => self.first.flags().is_mate_reverse_complemented(),
        }
    }

    /// The first read's distance from its 5′ end, or else the second read's template end
    /// distance.
    pub(crate) fn five_prime_distance(&self) -> Result<Option<usize>> {
        if let Some(distance) = self
            .first_read()
            .and_then(TemplateRead::five_prime_distance)
        {
            return Ok(Some(distance));
        }
        match self.second_read() {
            Some(second) => second.template_end_distance(),
            None => Ok(None),
        }
    }

    /// The first read's template end distance, or else the second read's distance from its 5′
    /// end where that read is of an FR pair.
    pub(crate) fn template_end_distance(&self) -> Result<Option<usize>> {
        match self.first_read() {
            Some(first) => first.template_end_distance(),
            None if self.first.is_fr_pair()? => Ok(self.first.five_prime_distance()),
            None => Ok(None),
        }
    }

    /// The template's first read here: the first of a pair, or a fragment's only read.
    fn first_read(&self) -> Option<&E> {
        self.entries().find(|entry| !is_second(entry.flags()))
    }

    /// The template's second read here: the last of a pair.
    fn second_read(&self) -> Option<&E> {
        self.entries().find(|entry| is_second(entry.flags()))
    }
}

/// The templates of a pileup's entries, each entry given with its template's number, from
/// [`number_templates`], and what it holds, which its template's vote counts.
pub(crate) fn templates<E>(
    entries: impl IntoIterator<Item = (usize, E, Observation)>,
    agreement: AgreementStrategy,
    disagreement: DisagreementStrategy,
    min_base_quality: i64,
) -> Vec<Template<E>> {
    let mut templates: Vec<Template<E>> = Vec::new();
    for (number, entry, observation) in entries {
        if let Some(template) = templates.get_mut(number) {
            template.push(entry);
        } else {
            templates.push(Template::of(entry));
        }
        templates[number]
            .vote
            .add(observation, agreement, disagreement, min_base_quality);
    }
    templates
}

/// One template at one pileup position: the reads of one name, their bases called into one.
///
/// A template's strand and distances are those of its first read, the first of a pair or a
/// fragment's only read, worked out from its second read where the first holds no base here.
#[derive(Debug)]
pub struct PileupTemplate<'a, R = bam::Record>(Template<PileupEntry<'a, R>>);

impl<R> Clone for PileupTemplate<'_, R> {
    fn clone(&self) -> Self {
        Self(self.0.clone())
    }
}

impl<'a, R: AlignmentRecord> PileupTemplate<'a, R> {
    /// The template's name, `*` for reads with none.
    pub fn name(&self) -> &'a BStr {
        name_of(self.0.first().record()).as_bstr()
    }

    /// The entries of the template's reads here, usually one or two, in input order.
    pub fn entries(&self) -> impl Iterator<Item = PileupEntry<'a, R>> + '_ {
        self.0.entries().copied()
    }

    /// What the template holds here: a base if a read's base votes, or else a deletion if a read's
    /// deletion votes, or else, with no vote, a base if any of its reads holds one, or else a
    /// deletion if any of them does, or else a skip.
    pub fn kind(&self) -> EntryKind {
        self.0.called().kind
    }

    /// Whether the template holds a deletion here, and no base that votes.
    pub fn is_deletion(&self) -> bool {
        self.kind() == EntryKind::Deletion
    }

    /// Whether every read of the template skips over the position.
    pub fn is_skip(&self) -> bool {
        self.kind() == EntryKind::Skip
    }

    /// Whether the template's base is a no-call, `N`.
    pub fn is_no_call(&self) -> bool {
        self.base() == Some(b'N')
    }

    /// The template's upper-cased base, its voting reads' bases called into one, or `None` without
    /// one.
    pub fn base(&self) -> Option<u8> {
        self.0.called().base
    }

    /// The quality of the template's base, or for a deletion of the next base, or `None` without a
    /// vote.
    pub fn quality(&self) -> Option<u8> {
        self.0.called().quality
    }

    /// Whether the template has a quality at the floor: a base, or a deletion followed by a base.
    pub fn passes(&self, min_base_quality: u8) -> bool {
        self.quality()
            .is_some_and(|quality| quality >= min_base_quality)
    }

    /// Whether the template's first read is aligned to the reverse strand: `false` for an F1R2
    /// pair and `true` for an F2R1 pair, read from the second read's flags without the first.
    pub fn is_reverse(&self) -> bool {
        self.0.is_reverse()
    }

    /// The number of the template's bases between its first read's 5′ end and this position: 0
    /// at that end.
    ///
    /// It is the first read's [`five_prime_distance`](PileupEntry::five_prime_distance) where it
    /// holds a base here, and otherwise the second read's
    /// [`template_end_distance`](PileupEntry::template_end_distance), which is an error for a read
    /// of an FR pair without a usable `MC` tag.
    pub fn five_prime_distance(&self) -> Result<Option<usize>> {
        self.0.five_prime_distance()
    }

    /// The number of the template's bases between this position and its other end, the 5′ end of
    /// the second read of an FR pair: 0 at that end.
    ///
    /// It is the first read's [`template_end_distance`](PileupEntry::template_end_distance), an
    /// error for a read of an FR pair without a usable `MC` tag, and without the first read here
    /// the second read's [`five_prime_distance`](PileupEntry::five_prime_distance) where that read
    /// [`is_fr_pair`](PileupEntry::is_fr_pair). It is `None` for any pair that is not FR.
    pub fn template_end_distance(&self) -> Result<Option<usize>> {
        self.0.template_end_distance()
    }
}

/// Whether a read is the last segment of a template.
fn is_second(flags: Flags) -> bool {
    flags.is_segmented() && flags.is_last_segment()
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

/// The name of a record, or `*` for a record with none.
pub(crate) fn record_name<R: noodles::sam::alignment::Record + ?Sized>(record: &R) -> String {
    record
        .name()
        .map_or_else(|| "*".to_owned(), |name| name.to_str_lossy().into_owned())
}
