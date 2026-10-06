use std::collections::{HashMap, VecDeque};
use std::io::{self, Read};
use std::mem;

use bstr::{BStr, ByteSlice};
use noodles::bam;
use noodles::sam::{self, alignment::record::Cigar as _, alignment::record::Flags};

use crate::auxiliary::{self, AuxValue};
use crate::error::{Error, Result};
use crate::footprint::Located;
use crate::pileup::{
    EntryKind, LiveRecord, NONE, OtherEnd, Pileup, RawEntry, quality_of, record_name,
};

/// Secondary, QC-fail, duplicate, and supplementary reads, which are left out by default.
pub const DEFAULT_EXCLUDE_FLAGS: Flags = Flags::from_bits_retain(0xF00);

/// The default minimum base quality of a pileup's filtered views, as in pysam's `pileup()` and
/// `samtools mpileup`.
pub const DEFAULT_MIN_BASE_QUALITY: u8 = 13;

const UNPLACED: usize = usize::MAX;

/// A coordinate-sorted stream of BAM records, which may be an unindexed pipe.
pub trait RecordSource {
    /// Reads the next record into `record`, reusing its buffer, and returns `false` at the end.
    fn read_record(&mut self, record: &mut bam::Record) -> io::Result<bool>;
}

impl<R: Read> RecordSource for bam::io::Reader<R> {
    fn read_record(&mut self, record: &mut bam::Record) -> io::Result<bool> {
        bam::io::Reader::read_record(self, record).map(|read| read > 0)
    }
}

impl<S: RecordSource + ?Sized> RecordSource for &mut S {
    fn read_record(&mut self, record: &mut bam::Record) -> io::Result<bool> {
        (**self).read_record(record)
    }
}

impl<S: RecordSource + ?Sized> RecordSource for Box<S> {
    fn read_record(&mut self, record: &mut bam::Record) -> io::Result<bool> {
        (**self).read_record(record)
    }
}

/// A [`RecordSource`] over any iterator of records.
#[derive(Debug)]
pub struct Records<I>(I);

impl<I> Records<I> {
    /// Wraps an iterator of records.
    pub fn new(records: I) -> Self {
        Self(records)
    }
}

impl<I: Iterator<Item = io::Result<bam::Record>>> RecordSource for Records<I> {
    fn read_record(&mut self, record: &mut bam::Record) -> io::Result<bool> {
        match self.0.next() {
            Some(next) => {
                *record = next?;
                Ok(true)
            }
            None => Ok(false),
        }
    }
}

type ReadFilter<'f> = Box<dyn FnMut(&bam::Record) -> bool + Send + 'f>;
type Tap<'f> = Box<dyn FnMut(bam::Record) -> io::Result<()> + Send + 'f>;

#[derive(Clone, Debug)]
struct Options {
    min_mapping_quality: u8,
    exclude_flags: Flags,
    min_base_quality: u8,
    proper_pairs_only: bool,
    without_overlaps: bool,
    aux_tags: Vec<[u8; 2]>,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            min_mapping_quality: 0,
            exclude_flags: DEFAULT_EXCLUDE_FLAGS,
            min_base_quality: DEFAULT_MIN_BASE_QUALITY,
            proper_pairs_only: false,
            without_overlaps: false,
            aux_tags: Vec::new(),
        }
    }
}

/// Whether a builder takes pileups, is handing the rest of its records to the tap, or is done.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Stage {
    Open,
    Closing,
    Closed,
}

#[derive(Clone, Copy, Debug, Default)]
struct TemplateMark {
    stamp: u32,
    first: u32,
    passing: u32,
}

/// Builds pileups from coordinate-sorted records in one forward pass.
///
/// Ask for pileups at positions that never move backwards: the same position again returns the
/// pileup already built, and an earlier one is an [`Error::Backwards`]. Positions are 0-based.
///
/// The builder holds a window of the records that span the position asked for, each with its
/// CIGAR decoded once into a footprint, and drops a record as soon as the builder has moved past
/// it. Pileups borrow the records in the window, and building one allocates nothing per entry.
///
/// Every record, filtered or not, is handed to a [`tap`](StreamingPileupBuilder::tap) exactly
/// once and in input order, as soon as the builder has moved past it and every record before it,
/// or when the builder closes. Keeping input order means a record is held until every record
/// before it has been passed, so a long record, such as one with a long reference skip, holds back
/// every record that starts within it. Without a tap, a record is dropped, and its buffers reused,
/// as soon as the builder has moved past it.
///
/// ```no_run
/// use noodles::bam;
/// use streampile::StreamingPileupBuilder;
///
/// let mut reader = bam::io::reader::Builder::default().build_from_path("reads.bam")?;
/// let header = reader.read_header()?;
/// let mut builder = StreamingPileupBuilder::new(reader, &header)?.min_mapping_quality(20);
///
/// let pileup = builder.pileup("chr1", 100)?;
/// println!("{} reads, {} at Q13 or more", pileup.unfiltered_depth(), pileup.filtered_depth());
///
/// let mut columns = builder.columns("chr1", 200, 300)?;
/// while let Some(pileup) = columns.next_pileup() {
///     let pileup = pileup?;
///     println!("{}\t{}", pileup.position(), pileup.bases().count());
/// }
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub struct StreamingPileupBuilder<'f, S: RecordSource> {
    source: S,
    header: sam::Header,
    options: Options,
    read_filter: Option<ReadFilter<'f>>,
    tap: Option<Tap<'f>>,
    pub(crate) slots: Vec<LiveRecord>,
    pub(crate) free: Vec<u32>,
    pub(crate) active: Vec<u32>,
    pub(crate) waiting: VecDeque<u32>,
    pub(crate) next: Option<u32>,
    exhausted: bool,
    active_reference_id: usize,
    min_active_end: i64,
    last_key: (usize, i64),
    at: Option<(usize, usize)>,
    built: bool,
    entries: Vec<RawEntry>,
    templates: HashMap<Vec<u8>, (u32, u32)>,
    free_templates: Vec<u32>,
    template_marks: Vec<TemplateMark>,
    stamp: u32,
    stage: Stage,
}

impl<'f, S: RecordSource> StreamingPileupBuilder<'f, S> {
    /// Starts a builder over the records of a source whose header declares coordinate order.
    ///
    /// No record is read until the first pileup is asked for.
    pub fn new(source: S, header: &sam::Header) -> Result<Self> {
        use sam::header::record::value::map::header::tag::SORT_ORDER;
        let sort_order = header
            .header()
            .and_then(|map| map.other_fields().get(&SORT_ORDER));
        if sort_order.is_none_or(|order| order.as_bytes() != b"coordinate") {
            return Err(Error::NotCoordinateSorted {
                found: sort_order.map(|order| order.to_str_lossy().into_owned()),
            });
        }
        Ok(Self {
            source,
            header: header.clone(),
            options: Options::default(),
            read_filter: None,
            tap: None,
            slots: Vec::new(),
            free: Vec::new(),
            active: Vec::new(),
            waiting: VecDeque::new(),
            next: None,
            exhausted: false,
            active_reference_id: UNPLACED,
            min_active_end: i64::MAX,
            last_key: (0, i64::MIN),
            at: None,
            built: false,
            entries: Vec::new(),
            templates: HashMap::new(),
            free_templates: Vec::new(),
            template_marks: Vec::new(),
            stamp: 0,
            stage: Stage::Open,
        })
    }

    /// Piles up only records with at least this mapping quality: 0 by default. A record with no
    /// mapping quality (255) passes every floor.
    #[must_use]
    pub fn min_mapping_quality(mut self, min_mapping_quality: u8) -> Self {
        self.options.min_mapping_quality = min_mapping_quality;
        self
    }

    /// Leaves out records with any of these flags: by default, secondary, QC-fail, duplicate,
    /// and supplementary records, [`DEFAULT_EXCLUDE_FLAGS`].
    #[must_use]
    pub fn exclude_flags(mut self, exclude_flags: Flags) -> Self {
        self.options.exclude_flags = exclude_flags;
        self
    }

    /// Sets the quality floor of each pileup's filtered views: [`DEFAULT_MIN_BASE_QUALITY`] by
    /// default. Entries under it are kept, and only the views leave them out.
    #[must_use]
    pub fn min_base_quality(mut self, min_base_quality: u8) -> Self {
        self.options.min_base_quality = min_base_quality;
        self
    }

    /// Piles up only records flagged as in a proper pair.
    #[must_use]
    pub fn proper_pairs_only(mut self, proper_pairs_only: bool) -> Self {
        self.options.proper_pairs_only = proper_pairs_only;
        self
    }

    /// Keeps one record per template, by name, in every pileup.
    ///
    /// The record kept is the first of its name whose entry is a base, or a deletion followed by
    /// a base, at the quality floor, or else the first of its name, and every entry of it is
    /// kept, its insertion entry included. Entries come in input order, so of two passing mates
    /// in a coordinate-sorted stream, the one that starts first is kept. A mate's skip, or its
    /// base under the floor, therefore never hides the other mate's base.
    #[must_use]
    pub fn without_overlaps(mut self, without_overlaps: bool) -> Self {
        self.options.without_overlaps = without_overlaps;
        self
    }

    /// Finds these auxiliary fields once, when each record is read, so entries reach them without
    /// a search, as for per-base tags read at every position of a deep pileup.
    #[must_use]
    pub fn index_aux_tags(mut self, tags: impl IntoIterator<Item = [u8; 2]>) -> Self {
        self.options.aux_tags = tags.into_iter().collect();
        self
    }

    /// Keeps a record for pileups only when this returns `true`. It is asked only of records that
    /// pass the other filters and align to the contig piled up, and a record it rejects still
    /// goes to the tap.
    #[must_use]
    pub fn read_filter(
        mut self,
        read_filter: impl FnMut(&bam::Record) -> bool + Send + 'f,
    ) -> Self {
        self.read_filter = Some(Box::new(read_filter));
        self
    }

    /// Hands every record to `tap` once the builder has moved past it, in input order.
    ///
    /// With a tap, [`close`](StreamingPileupBuilder::close) reads the rest of the input, so an
    /// output written by the tap is complete.
    #[must_use]
    pub fn tap(mut self, tap: impl FnMut(bam::Record) -> io::Result<()> + Send + 'f) -> Self {
        self.tap = Some(Box::new(tap));
        self
    }

    /// The header the records were read with.
    pub fn header(&self) -> &sam::Header {
        &self.header
    }

    /// Advances to a position on a contig, at or after the last one, and piles up the records
    /// there.
    pub fn pileup(&mut self, contig: &str, position: usize) -> Result<Pileup<'_>> {
        if self.stage != Stage::Open {
            return Err(Error::Closed);
        }
        let reference_sequence_id = self.reference_sequence_id(contig)?;
        self.pileup_at(reference_sequence_id, position)
    }

    /// Advances to a position on the contig with this index in the header, at or after the last
    /// one, and piles up the records there.
    pub fn pileup_at(
        &mut self,
        reference_sequence_id: usize,
        position: usize,
    ) -> Result<Pileup<'_>> {
        if self.stage != Stage::Open {
            return Err(Error::Closed);
        }
        if self.built && self.at == Some((reference_sequence_id, position)) {
            return Ok(self.view(reference_sequence_id, position));
        }
        if reference_sequence_id >= self.header.reference_sequences().len() {
            return Err(Error::UnknownReferenceSequenceId(reference_sequence_id));
        }
        if let Some((from_id, from_position)) = self.at
            && (reference_sequence_id, position) < (from_id, from_position)
        {
            return Err(Error::Backwards {
                contig: self.name_of(reference_sequence_id).to_string(),
                position,
                from_contig: self.name_of(from_id).to_string(),
                from_position,
            });
        }
        self.at = Some((reference_sequence_id, position));
        self.built = false;
        self.advance(reference_sequence_id, position as i64)?;
        self.collect(position as i64);
        self.built = true;
        Ok(self.view(reference_sequence_id, position))
    }

    /// Sweeps the pileup at every position of a span of a contig, covered or not.
    pub fn columns(
        &mut self,
        contig: &str,
        start: usize,
        end: usize,
    ) -> Result<Columns<'_, 'f, S>> {
        let reference_sequence_id = self.reference_sequence_id(contig)?;
        self.columns_at(reference_sequence_id, start, end)
    }

    /// Sweeps the pileup at every position of a span of the contig with this index in the
    /// header, covered or not.
    pub fn columns_at(
        &mut self,
        reference_sequence_id: usize,
        start: usize,
        end: usize,
    ) -> Result<Columns<'_, 'f, S>> {
        if end < start {
            return Err(Error::InvalidSpan { start, end });
        }
        Ok(Columns {
            builder: self,
            reference_sequence_id,
            next: start,
            end,
        })
    }

    /// Stops, first handing every record not yet handed to the tap to it, in input order.
    ///
    /// With a tap, the rest of the input is read to the end, so an output written by the tap is
    /// complete. Without one, no more of the input is read. A record goes to the tap once even
    /// when the tap fails, and closing again after an error hands over the records after it.
    /// Dropping a builder closes it too, but ignores any error, so close it to see one.
    pub fn close(&mut self) -> Result<()> {
        if self.stage == Stage::Closed {
            return Ok(());
        }
        self.stage = Stage::Closing;
        self.built = false;
        self.active.clear();
        self.entries.clear();
        if let Some(tap) = self.tap.as_mut() {
            while let Some(index) = self.waiting.pop_front() {
                tap(mem::take(&mut self.slots[index as usize].record))?;
            }
            if let Some(index) = self.next.take() {
                tap(mem::take(&mut self.slots[index as usize].record))?;
            }
            if !self.exhausted {
                let mut record = bam::Record::default();
                while self.source.read_record(&mut record)? {
                    tap(mem::take(&mut record))?;
                }
                self.exhausted = true;
            }
        }
        self.slots.clear();
        self.free.clear();
        self.waiting.clear();
        self.next = None;
        self.stage = Stage::Closed;
        Ok(())
    }

    fn reference_sequence_id(&self, contig: &str) -> Result<usize> {
        self.header
            .reference_sequences()
            .get_index_of(contig.as_bytes().as_bstr())
            .ok_or_else(|| Error::UnknownContig(contig.to_owned()))
    }

    fn name_of(&self, reference_sequence_id: usize) -> &BStr {
        self.header
            .reference_sequences()
            .get_index(reference_sequence_id)
            .map_or(b"".as_bstr(), |(name, _)| name.as_bstr())
    }

    fn view(&self, reference_sequence_id: usize, position: usize) -> Pileup<'_> {
        Pileup {
            reference_sequence_id,
            reference_sequence_name: self.name_of(reference_sequence_id),
            position,
            min_base_quality: self.options.min_base_quality,
            entries: &self.entries,
            slots: &self.slots,
            aux_tags: &self.options.aux_tags,
        }
    }

    fn advance(&mut self, reference_sequence_id: usize, pos: i64) -> Result<()> {
        if reference_sequence_id != self.active_reference_id {
            self.evict(None);
            self.active_reference_id = reference_sequence_id;
        } else if self.min_active_end <= pos {
            self.evict(Some(pos));
        }
        let target = (reference_sequence_id, pos + 1);
        while let Some(index) = self.peek()? {
            let live = &self.slots[index as usize];
            let key = (live.reference_id, live.start);
            if key > target {
                break;
            }
            if key < self.last_key {
                return Err(Error::OutOfOrder {
                    name: record_name(&live.record),
                });
            }
            self.last_key = key;
            self.next = None;
            let accepted = if key.0 == reference_sequence_id {
                self.accept(index, pos)
            } else {
                Ok(false)
            };
            if self.tap.is_some() {
                self.waiting.push_back(index);
            }
            if let Ok(true) = accepted {
                self.activate(index);
            } else if self.tap.is_none() {
                self.free.push(index);
            }
            accepted?;
        }
        self.release(reference_sequence_id, pos)
    }

    fn peek(&mut self) -> Result<Option<u32>> {
        if self.next.is_some() || self.exhausted {
            return Ok(self.next);
        }
        let index = if let Some(index) = self.free.pop() {
            index
        } else {
            self.slots.push(LiveRecord::default());
            u32::try_from(self.slots.len() - 1)
                .map_err(|_| io::Error::other("more than 4,294,967,295 records are held"))?
        };
        let tapped = self.tap.is_some();
        let live = &mut self.slots[index as usize];
        match self.source.read_record(&mut live.record) {
            Ok(true) => {
                if let Err(error) = describe(live, tapped) {
                    live.reference_id = UNPLACED;
                    if tapped {
                        self.waiting.push_back(index);
                    } else {
                        self.free.push(index);
                    }
                    return Err(error.into());
                }
                self.next = Some(index);
            }
            read => {
                self.free.push(index);
                read?;
                self.exhausted = true;
            }
        }
        Ok(self.next)
    }

    fn accept(&mut self, index: u32, pos: i64) -> Result<bool> {
        let live = &mut self.slots[index as usize];
        let flags = live.flags;
        let mapping_quality = live.record.mapping_quality().map_or(255, u8::from);
        if flags.intersects(self.options.exclude_flags)
            || mapping_quality < self.options.min_mapping_quality
            || (self.options.proper_pairs_only && !flags.is_properly_segmented())
            || flags.is_unmapped()
            || live.start < 0
        {
            return Ok(false);
        }
        let stored_bases = live.record.sequence().len();
        let placed = live
            .footprint
            .fill(live.start, live.record.cigar().as_bytes(), stored_bases)
            .map_err(|source| Error::InvalidRecord {
                name: record_name(&live.record),
                source,
            })?;
        if !placed {
            return Ok(false);
        }
        if let Some(read_filter) = self.read_filter.as_mut()
            && !read_filter(&live.record)
        {
            return Ok(false);
        }
        if live.footprint.end <= pos {
            return Ok(false);
        }
        index_fields(live, &self.options.aux_tags).map_err(|error| match error {
            Error::Io(source) => Error::InvalidRecord {
                name: record_name(&live.record),
                source,
            },
            error => error,
        })?;
        Ok(true)
    }

    fn activate(&mut self, index: u32) {
        let end = self.slots[index as usize].footprint.end;
        self.min_active_end = self.min_active_end.min(end);
        self.active.push(index);
        if self.options.without_overlaps {
            self.register(index);
        }
    }

    fn evict(&mut self, before: Option<i64>) {
        let mut active = mem::take(&mut self.active);
        let mut kept = 0;
        let mut min_end = i64::MAX;
        for at in 0..active.len() {
            let index = active[at];
            let end = self.slots[index as usize].footprint.end;
            if before.is_some_and(|pos| end > pos) {
                active[kept] = index;
                kept += 1;
                min_end = min_end.min(end);
            } else {
                self.unregister(index);
                if self.tap.is_none() {
                    self.free.push(index);
                }
            }
        }
        active.truncate(kept);
        self.active = active;
        self.min_active_end = min_end;
    }

    fn release(&mut self, reference_sequence_id: usize, pos: i64) -> Result<()> {
        let Some(tap) = self.tap.as_mut() else {
            return Ok(());
        };
        while let Some(&index) = self.waiting.front() {
            let live = &mut self.slots[index as usize];
            let passed = live.reference_id < reference_sequence_id
                || (live.reference_id == reference_sequence_id && live.end <= pos);
            if !passed {
                break;
            }
            self.waiting.pop_front();
            let record = mem::take(&mut live.record);
            self.free.push(index);
            tap(record)?;
        }
        Ok(())
    }

    fn collect(&mut self, pos: i64) {
        self.entries.clear();
        for &index in &self.active {
            let footprint = &mut self.slots[index as usize].footprint;
            let entry = |kind, offset| RawEntry {
                slot: index,
                kind,
                offset,
                length: 0,
            };
            match footprint.locate(pos) {
                Some(Located::Base(offset)) => self.entries.push(entry(EntryKind::Base, offset)),
                Some(Located::Deletion(next)) => {
                    self.entries
                        .push(entry(EntryKind::Deletion, next.unwrap_or(NONE)));
                }
                Some(Located::Skip) => self.entries.push(entry(EntryKind::Skip, NONE)),
                None => {}
            }
            if let Some(insertion) = footprint.insertion_at(pos) {
                self.entries.push(RawEntry {
                    slot: index,
                    kind: EntryKind::Insertion,
                    offset: insertion.offset,
                    length: insertion.length,
                });
            }
        }
        if self.options.without_overlaps {
            self.keep_one_per_template();
        }
    }

    fn keep_one_per_template(&mut self) {
        self.stamp = self.stamp.wrapping_add(1);
        if self.stamp == 0 {
            self.template_marks.fill(TemplateMark::default());
            self.stamp = 1;
        }
        let stamp = self.stamp;
        let floor = self.options.min_base_quality;
        for raw in &self.entries {
            let live = &self.slots[raw.slot as usize];
            let mark = &mut self.template_marks[live.template as usize];
            if mark.stamp != stamp {
                *mark = TemplateMark {
                    stamp,
                    first: raw.slot,
                    passing: NONE,
                };
            }
            if mark.passing == NONE
                && quality_of(&live.record, raw.kind, raw.offset).is_some_and(|q| q >= floor)
            {
                mark.passing = raw.slot;
            }
        }
        let (slots, marks) = (&self.slots, &self.template_marks);
        self.entries.retain(|raw| {
            let mark = marks[slots[raw.slot as usize].template as usize];
            raw.slot
                == if mark.passing == NONE {
                    mark.first
                } else {
                    mark.passing
                }
        });
    }

    fn register(&mut self, index: u32) {
        let live = &self.slots[index as usize];
        let name = live
            .record
            .name()
            .map_or(b"*".as_slice(), |name| name.as_bytes());
        let id = if let Some((id, count)) = self.templates.get_mut(name) {
            *count += 1;
            *id
        } else {
            let id = self.free_templates.pop().unwrap_or_else(|| {
                self.template_marks.push(TemplateMark::default());
                (self.template_marks.len() - 1) as u32
            });
            self.templates.insert(name.to_vec(), (id, 1));
            id
        };
        self.slots[index as usize].template = id;
    }

    fn unregister(&mut self, index: u32) {
        let live = &mut self.slots[index as usize];
        if live.template == NONE {
            return;
        }
        live.template = NONE;
        let name = live
            .record
            .name()
            .map_or(b"*".as_slice(), |name| name.as_bytes());
        if let Some((id, count)) = self.templates.get_mut(name) {
            *count -= 1;
            if *count == 0 {
                self.free_templates.push(*id);
                self.templates.remove(name);
            }
        }
    }
}

impl<S: RecordSource> Drop for StreamingPileupBuilder<'_, S> {
    fn drop(&mut self) {
        if !std::thread::panicking() {
            let _ = self.close();
        }
    }
}

/// A forward sweep of the pileup at every position of a span, from
/// [`columns`](StreamingPileupBuilder::columns).
pub struct Columns<'b, 'f, S: RecordSource> {
    builder: &'b mut StreamingPileupBuilder<'f, S>,
    reference_sequence_id: usize,
    next: usize,
    end: usize,
}

impl<S: RecordSource> Columns<'_, '_, S> {
    /// The pileup at the next position of the span, or `None` past its end.
    pub fn next_pileup(&mut self) -> Option<Result<Pileup<'_>>> {
        if self.next >= self.end {
            return None;
        }
        let position = self.next;
        self.next += 1;
        Some(self.builder.pileup_at(self.reference_sequence_id, position))
    }
}

/// Reads a record's contig, start, and flags, and with a tap, where the builder is past it.
fn describe(live: &mut LiveRecord, tapped: bool) -> io::Result<()> {
    let record = &live.record;
    live.reference_id = record
        .reference_sequence_id()
        .transpose()?
        .unwrap_or(UNPLACED);
    live.start = record
        .alignment_start()
        .transpose()?
        .map_or(-1, |start| usize::from(start) as i64 - 1);
    live.flags = record.flags();
    live.template = NONE;
    live.other_end = OtherEnd::None;
    if tapped {
        let cigar = record.cigar();
        let length = if live.flags.is_unmapped() || cigar.is_empty() {
            1
        } else {
            cigar.alignment_span()?.max(1) as i64
        };
        live.end = live.start + length;
    }
    Ok(())
}

/// Finds the indexed auxiliary fields of an accepted record, and the 5′ end of its mate.
fn index_fields(live: &mut LiveRecord, tags: &[[u8; 2]]) -> Result<()> {
    live.fields.clear();
    live.fields.resize(tags.len(), None);
    let mut mate_cigar = MateCigar::Unsearched;
    if !tags.is_empty() {
        let fields = &mut live.fields;
        let mut found = None;
        auxiliary::walk(live.record.data().as_bytes(), |tag, field| {
            if tag == *b"MC" {
                found.get_or_insert(field);
            }
            if let Some(index) = tags.iter().position(|wanted| *wanted == tag)
                && fields[index].is_none()
            {
                fields[index] = Some(field);
            }
        })?;
        mate_cigar = found.map_or(MateCigar::Missing, MateCigar::Found);
    }
    live.other_end = other_end(live, mate_cigar)?;
    Ok(())
}

/// Where a record's `MC` field is, if it has been looked for.
#[derive(Clone, Copy)]
enum MateCigar {
    Unsearched,
    Missing,
    Found(auxiliary::Field),
}

/// Where the 5′ end of the mate of a read in an FR pair is.
///
/// For a reverse read it is the mate's start. For a forward read it is the end of the mate's
/// alignment, from its start and its `MC` tag alone, never from the template length (TLEN).
fn other_end(live: &LiveRecord, mate_cigar: MateCigar) -> Result<OtherEnd> {
    let flags = live.flags;
    let reverse = flags.is_reverse_complemented();
    if !flags.is_segmented()
        || flags.is_mate_unmapped()
        || reverse == flags.is_mate_reverse_complemented()
    {
        return Ok(OtherEnd::None);
    }
    let record = &live.record;
    if record.mate_reference_sequence_id().transpose()? != Some(live.reference_id) {
        return Ok(OtherEnd::None);
    }
    let Some(mate_start) = record.mate_alignment_start().transpose()? else {
        return Ok(OtherEnd::None);
    };
    let mate_start = usize::from(mate_start) as i64 - 1;
    if reverse {
        return Ok(if mate_start < live.footprint.end {
            OtherEnd::At(mate_start)
        } else {
            OtherEnd::None
        });
    }
    let data = record.data().as_bytes();
    let mate_cigar = match mate_cigar {
        MateCigar::Found(field) => Some(field.value(data)?),
        MateCigar::Missing => None,
        MateCigar::Unsearched => auxiliary::find(data, *b"MC")?,
    };
    Ok(match mate_cigar {
        None => OtherEnd::MissingMateCigar,
        Some(value) => match mate_span(value) {
            None => OtherEnd::InvalidMateCigar,
            Some(span) if live.start < mate_start + span => OtherEnd::At(mate_start + span - 1),
            Some(_) => OtherEnd::None,
        },
    })
}

/// The number of reference bases an `MC` value spans, or `None` for a value that is not a
/// CIGAR string spanning at least one base with operators no longer than BAM allows.
pub(crate) fn mate_span(value: AuxValue<'_>) -> Option<i64> {
    const LONGEST_OPERATOR: usize = (1 << 28) - 1;
    let AuxValue::String(text) = value else {
        return None;
    };
    let cigar = sam::record::Cigar::new(text);
    if cigar
        .iter()
        .any(|op| op.map_or(true, |op| op.len() > LONGEST_OPERATOR))
    {
        return None;
    }
    let span = cigar.alignment_span().ok()?;
    i64::try_from(span).ok().filter(|&span| span > 0)
}
