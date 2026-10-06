use std::collections::VecDeque;
use std::io;
use std::mem;

use bstr::{BStr, ByteSlice};
use noodles::bam;
use noodles::sam::alignment::record::cigar::Op;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::{self, alignment::record::Cigar as _, alignment::record::Flags};

use crate::auxiliary::{self, AuxValue};
use crate::error::{Error, Result};
use crate::footprint::{Footprint, Located};
use crate::pileup::{
    EntryKind, LiveRecord, NONE, OtherEnd, Pileup, RawEntry, name_hash, quality_of, record_name,
    templates_kept,
};
use crate::source::{AlignmentRecord, RecordSource};

/// Secondary, QC-fail, duplicate, and supplementary reads, which are left out by default.
pub const DEFAULT_EXCLUDE_FLAGS: Flags = Flags::from_bits_retain(0xF00);

/// The default minimum base quality of a pileup's filtered views, as in pysam's `pileup()` and
/// `samtools mpileup`.
pub const DEFAULT_MIN_BASE_QUALITY: u8 = 13;

const UNPLACED: usize = usize::MAX;

type ReadFilter<'f, R> = Box<dyn FnMut(&R) -> io::Result<bool> + Send + 'f>;
type Tap<'f, R> = Box<dyn FnMut(R) -> io::Result<()> + Send + 'f>;

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
    read_filter: Option<ReadFilter<'f, S::Record>>,
    tap: Option<Tap<'f, S::Record>>,
    pub(crate) slots: Vec<LiveRecord<S::Record>>,
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
    pub fn read_filter(self, mut read_filter: impl FnMut(&S::Record) -> bool + Send + 'f) -> Self {
        self.try_read_filter(move |record| Ok(read_filter(record)))
    }

    /// Keeps a record for pileups only when this returns `Ok(true)`, as
    /// [`read_filter`](StreamingPileupBuilder::read_filter) does, and stops at an error.
    #[must_use]
    pub(crate) fn try_read_filter(
        mut self,
        read_filter: impl FnMut(&S::Record) -> io::Result<bool> + Send + 'f,
    ) -> Self {
        self.read_filter = Some(Box::new(read_filter));
        self
    }

    /// Hands every record to `tap` once the builder has moved past it, in input order.
    ///
    /// With a tap, [`close`](StreamingPileupBuilder::close) reads the rest of the input, so an
    /// output written by the tap is complete.
    #[must_use]
    pub fn tap(mut self, tap: impl FnMut(S::Record) -> io::Result<()> + Send + 'f) -> Self {
        self.tap = Some(Box::new(tap));
        self
    }

    /// The header the records were read with.
    pub fn header(&self) -> &sam::Header {
        &self.header
    }

    /// Whether a record passes the flag, mapping-quality, and proper-pair filters, is placed
    /// with a reference-consuming CIGAR operator, and then passes the read filter, so that it is
    /// piled up wherever it spans.
    pub fn accepts(&mut self, record: &S::Record) -> Result<bool> {
        let bam = record.bam();
        let invalid = |source| Error::InvalidRecord {
            name: record_name(bam),
            source,
        };
        let start = bam.alignment_start().transpose().map_err(invalid)?;
        let placed = bam.reference_sequence_id().is_some()
            && self.passes(bam.flags(), bam.mapping_quality().map_or(255, u8::from))
            && Footprint::default()
                .fill(
                    start.map_or(-1, |start| usize::from(start) as i64 - 1),
                    bam.cigar().as_bytes(),
                    bam.sequence().len(),
                )
                .map_err(invalid)?;
        Ok(placed
            && start.is_some()
            && match self.read_filter.as_mut() {
                Some(read_filter) => read_filter(record)?,
                None => true,
            })
    }

    /// Advances to a position on a contig, at or after the last one, and piles up the records
    /// there.
    pub fn pileup(&mut self, contig: &str, position: usize) -> Result<Pileup<'_, S::Record>> {
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
    ) -> Result<Pileup<'_, S::Record>> {
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
                let mut record = S::Record::default();
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

    /// Stops without handing any more records to the tap, so dropping the builder runs no tap.
    #[cfg(feature = "python")]
    pub(crate) fn abandon(&mut self) {
        self.stage = Stage::Closed;
        self.built = false;
        self.active.clear();
        self.entries.clear();
        self.slots.clear();
        self.free.clear();
        self.waiting.clear();
        self.next = None;
    }

    fn passes(&self, flags: Flags, mapping_quality: u8) -> bool {
        !(flags.intersects(self.options.exclude_flags)
            || mapping_quality < self.options.min_mapping_quality
            || (self.options.proper_pairs_only && !flags.is_properly_segmented())
            || flags.is_unmapped())
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

    fn view(&self, reference_sequence_id: usize, position: usize) -> Pileup<'_, S::Record> {
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

    fn free(&mut self, index: u32) {
        self.slots[index as usize].record.release();
        self.free.push(index);
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
                    name: record_name(live.record.bam()),
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
                self.free(index);
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
                        self.free(index);
                    }
                    return Err(error.into());
                }
                self.next = Some(index);
            }
            read => {
                self.free(index);
                read?;
                self.exhausted = true;
            }
        }
        Ok(self.next)
    }

    fn accept(&mut self, index: u32, pos: i64) -> Result<bool> {
        let live = &self.slots[index as usize];
        let mapping_quality = live.record.bam().mapping_quality().map_or(255, u8::from);
        if !self.passes(live.flags, mapping_quality) || live.start < 0 {
            return Ok(false);
        }
        let live = &mut self.slots[index as usize];
        let record = live.record.bam();
        let invalid = |source| Error::InvalidRecord {
            name: record_name(record),
            source,
        };
        let placed = live
            .footprint
            .fill(
                live.start,
                record.cigar().as_bytes(),
                record.sequence().len(),
            )
            .map_err(invalid)?;
        if !placed {
            return Ok(false);
        }
        if let Some(read_filter) = self.read_filter.as_mut()
            && !read_filter(&live.record)?
        {
            return Ok(false);
        }
        if live.footprint.end <= pos {
            return Ok(false);
        }
        live.name_hash = name_hash(live.record.bam());
        index_fields(live, &self.options.aux_tags).map_err(|error| match error {
            Error::Io(source) => Error::InvalidRecord {
                name: record_name(live.record.bam()),
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
            } else if self.tap.is_none() {
                self.free(index);
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
        let floor = self.options.min_base_quality;
        let slots = &self.slots;
        let kept = templates_kept(&self.entries, |raw| {
            let live = &slots[raw.slot as usize];
            let quality = quality_of(live.record.bam(), raw.kind, raw.offset);
            (
                live.template(),
                raw.slot as usize,
                quality.is_some_and(|q| q >= floor),
            )
        });
        let mut kept = kept.into_iter();
        self.entries.retain(|_| kept.next().unwrap_or(true));
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
    pub fn next_pileup(&mut self) -> Option<Result<Pileup<'_, S::Record>>> {
        if self.next >= self.end {
            return None;
        }
        let position = self.next;
        self.next += 1;
        Some(self.builder.pileup_at(self.reference_sequence_id, position))
    }
}

/// Reads a record's contig, start, and flags, and with a tap, where the builder is past it.
fn describe<R: AlignmentRecord>(live: &mut LiveRecord<R>, tapped: bool) -> io::Result<()> {
    let record = live.record.bam();
    live.reference_id = record
        .reference_sequence_id()
        .transpose()?
        .unwrap_or(UNPLACED);
    live.start = record
        .alignment_start()
        .transpose()?
        .map_or(-1, |start| usize::from(start) as i64 - 1);
    live.flags = record.flags();
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
fn index_fields<R: AlignmentRecord>(live: &mut LiveRecord<R>, tags: &[[u8; 2]]) -> Result<()> {
    live.fields.clear();
    live.fields.resize(tags.len(), None);
    let record = live.record.bam();
    let mut mate_cigar = MateCigar::Unsearched;
    if !tags.is_empty() {
        let fields = &mut live.fields;
        let mut found = None;
        auxiliary::walk(record.data().as_bytes(), |tag, field| {
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
    live.other_end = other_end(
        record,
        live.reference_id,
        live.start,
        live.footprint.end,
        mate_cigar,
    )?;
    Ok(())
}

/// Where a record's `MC` field is, if it has been looked for.
#[derive(Clone, Copy)]
pub(crate) enum MateCigar {
    Unsearched,
    Missing,
    Found(auxiliary::Field),
}

/// Where the unclipped 5′ end of the mate of a read in an FR pair is, for a read placed from
/// `start` to `end` on the contig with index `reference_id`.
///
/// For a reverse read it is the mate's start less the clips before it, and for a forward read the
/// end of the mate's alignment plus the clips after it, both from the mate's start and its `MC`
/// tag alone, never from the template length (TLEN).
pub(crate) fn other_end(
    record: &bam::Record,
    reference_id: usize,
    start: i64,
    end: i64,
    mate_cigar: MateCigar,
) -> Result<OtherEnd> {
    let flags = record.flags();
    let reverse = flags.is_reverse_complemented();
    if !flags.is_segmented()
        || flags.is_mate_unmapped()
        || reverse == flags.is_mate_reverse_complemented()
    {
        return Ok(OtherEnd::None);
    }
    if record.mate_reference_sequence_id().transpose()? != Some(reference_id) {
        return Ok(OtherEnd::None);
    }
    let Some(mate_start) = record.mate_alignment_start().transpose()? else {
        return Ok(OtherEnd::None);
    };
    let mate_start = usize::from(mate_start) as i64 - 1;
    let data = record.data().as_bytes();
    let mate_cigar = match mate_cigar {
        MateCigar::Found(field) => Some(field.value(data)?),
        MateCigar::Missing => None,
        MateCigar::Unsearched => auxiliary::find(data, *b"MC")?,
    };
    let Some(value) = mate_cigar else {
        return Ok(OtherEnd::MissingMateCigar);
    };
    let Some(extent) = mate_extent(value) else {
        return Ok(OtherEnd::InvalidMateCigar);
    };
    let five_prime = if reverse {
        mate_start - extent.leading
    } else {
        mate_start + extent.span - 1 + extent.trailing
    };
    let reachable = if reverse {
        five_prime < end
    } else {
        start <= five_prime
    };
    Ok(if reachable {
        OtherEnd::At(five_prime)
    } else {
        OtherEnd::None
    })
}

/// The extent of a mate's alignment from its `MC` value, in reference bases.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct MateExtent {
    /// The soft- and hard-clipped bases before the alignment.
    pub(crate) leading: i64,
    /// The reference bases the alignment spans.
    pub(crate) span: i64,
    /// The soft- and hard-clipped bases after the alignment.
    pub(crate) trailing: i64,
}

/// The extent of the alignment an `MC` value describes, or `None` for a value that is not a CIGAR
/// string spanning at least one base with operators no longer than BAM allows.
pub(crate) fn mate_extent(value: AuxValue<'_>) -> Option<MateExtent> {
    const LONGEST_OPERATOR: usize = (1 << 28) - 1;
    let AuxValue::String(text) = value else {
        return None;
    };
    let cigar = sam::record::Cigar::new(text);
    let ops = cigar
        .iter()
        .collect::<std::result::Result<Vec<_>, _>>()
        .ok()?;
    if ops.iter().any(|op| op.len() > LONGEST_OPERATOR) {
        return None;
    }
    let span = i64::try_from(cigar.alignment_span().ok()?)
        .ok()
        .filter(|&span| span > 0)?;
    let is_clip = |op: &&Op| matches!(op.kind(), Kind::SoftClip | Kind::HardClip);
    let clipped = |ops: &mut dyn Iterator<Item = &Op>| -> i64 {
        ops.take_while(is_clip).map(|op| op.len() as i64).sum()
    };
    Some(MateExtent {
        leading: clipped(&mut ops.iter()),
        span,
        trailing: clipped(&mut ops.iter().rev()),
    })
}
