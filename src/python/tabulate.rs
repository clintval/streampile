//! The per-read and per-base work of `tabulate`: alleles, runs, normalization, and ledgers.

use std::collections::{BTreeMap, HashMap};

use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use pyo3::exceptions::{PyIndexError, PyValueError};
use pyo3::intern;
use pyo3::prelude::*;
use pyo3::sync::PyOnceLock;
use pyo3::types::{PyDict, PyString, PyTuple};

use super::{Int, bridge};
use crate::DEFAULT_EXCLUDE_FLAGS;

const REFERENCE_PADDING: i64 = 10_000;

/// An allele's reference and alternate bases.
type AlleleKey = (Vec<u8>, Vec<u8>);

/// The reads of an allele on the forward and reverse strands.
type Strands = [i64; 2];

/// A normalized allele: its position and its reference and alternate bases.
type Normalized<T> = (i64, Vec<T>, Vec<T>);

static TABULATED_BASE: PyOnceLock<Py<PyAny>> = PyOnceLock::new();

/// Two sliding windows over one contig of an indexed FASTA, upper-cased.
///
/// A read with a long reference skip reads two far-apart stretches, so the window used last is
/// kept beside the other, and neither is fetched again for the next read.
struct Reference {
    fasta: Py<PyAny>,
    contig: Py<PyString>,
    length: i64,
    windows: Vec<(i64, i64, Vec<u8>)>,
}

impl Reference {
    fn new(fasta: &Bound<'_, PyAny>, contig: &Bound<'_, PyString>) -> PyResult<Self> {
        Ok(Self {
            fasta: fasta.clone().unbind(),
            contig: contig.clone().unbind(),
            length: fasta
                .call_method1("get_reference_length", (contig,))?
                .extract()?,
            windows: Vec::new(),
        })
    }

    /// The upper-cased bases from `start` to `end`, cut to the contig.
    fn get(&mut self, py: Python<'_>, start: i64, end: i64) -> PyResult<&[u8]> {
        let (start, end) = (start.max(0), end.min(self.length));
        let found = self
            .windows
            .iter()
            .position(|(low, high, _)| *low <= start && end <= *high);
        let index = if let Some(index) = found {
            if index > 0 {
                self.windows.reverse();
            }
            0
        } else {
            let low = (start - REFERENCE_PADDING).max(0);
            let high = (end.max(start) + 10 * REFERENCE_PADDING).min(self.length);
            let fetched = self
                .fasta
                .bind(py)
                .call_method1("fetch", (self.contig.bind(py), low, high))?;
            let mut bases = fetched.cast::<PyString>()?.to_str()?.as_bytes().to_vec();
            bases.make_ascii_uppercase();
            self.windows.truncate(1);
            self.windows.insert(0, (low, high, bases));
            0
        };
        let (low, _, bases) = &self.windows[index];
        Ok(slice(bases, start - low, end - low))
    }
}

/// A Python slice of bytes from `start` to `end`, with non-negative bounds cut to the bytes.
fn slice(bytes: &[u8], start: i64, end: i64) -> &[u8] {
    let length = bytes.len() as i64;
    let start = start.clamp(0, length) as usize;
    let end = end.clamp(0, length) as usize;
    &bytes[start..end.max(start)]
}

fn is_acgt(base: u8) -> bool {
    matches!(base, b'A' | b'C' | b'G' | b'T')
}

/// The counts of one territory span while reads are added to it, by strand.
struct Ledger {
    start: i64,
    end: i64,
    depth: Vec<i64>,
    ref_fwd: Vec<i64>,
    ref_rev: Vec<i64>,
    no_calls: Vec<i64>,
    alleles: HashMap<i64, HashMap<AlleleKey, Strands>>,
}

impl Ledger {
    fn new(start: i64, end: i64) -> Self {
        let length = (end - start).max(0) as usize;
        Self {
            start,
            end,
            depth: vec![0; length + 1],
            ref_fwd: vec![0; length + 1],
            ref_rev: vec![0; length + 1],
            no_calls: vec![0; length],
            alleles: HashMap::new(),
        }
    }

    fn add(&mut self, start: i64, end: i64, reference: bool, reverse: bool) {
        let (start, end) = (start.max(self.start), end.min(self.end));
        if start < end {
            let (from, to) = ((start - self.start) as usize, (end - self.start) as usize);
            self.depth[from] += 1;
            self.depth[to] -= 1;
            if reference {
                let counts = if reverse {
                    &mut self.ref_rev
                } else {
                    &mut self.ref_fwd
                };
                counts[from] += 1;
                counts[to] -= 1;
            }
        }
    }

    fn add_no_call(&mut self, pos: i64) {
        if self.start <= pos && pos < self.end {
            self.no_calls[(pos - self.start) as usize] += 1;
        }
    }

    fn add_allele(&mut self, allele: &Allele, reverse: bool) {
        if self.start <= allele.pos && allele.pos < self.end {
            let counts = self
                .alleles
                .entry(allele.pos)
                .or_default()
                .entry((allele.reference.clone(), allele.alternate.clone()))
                .or_default();
            counts[usize::from(reverse)] += 1;
        }
    }
}

/// The rows of a finished ledger, made one at a time as they are asked for.
struct Rows {
    ledger: Ledger,
    bases: Vec<u8>,
    offset: usize,
    depth: i64,
    ref_fwd: i64,
    ref_rev: i64,
}

impl Rows {
    fn new(py: Python<'_>, ledger: Ledger, reference: &mut Reference) -> PyResult<Self> {
        let bases = reference.get(py, ledger.start, ledger.end)?.to_vec();
        Ok(Self {
            ledger,
            bases,
            offset: 0,
            depth: 0,
            ref_fwd: 0,
            ref_rev: 0,
        })
    }

    fn next(
        &mut self,
        py: Python<'_>,
        contig: &Bound<'_, PyString>,
    ) -> PyResult<Option<Py<PyAny>>> {
        let Some(&base) = self.bases.get(self.offset) else {
            return Ok(None);
        };
        let offset = self.offset;
        self.offset += 1;
        let ledger = &mut self.ledger;
        self.depth += ledger.depth[offset];
        self.ref_fwd += ledger.ref_fwd[offset];
        self.ref_rev += ledger.ref_rev[offset];
        let pos = ledger.start + offset as i64;
        let fields = PyDict::new(py);
        fields.set_item(intern!(py, "contig"), contig)?;
        fields.set_item(intern!(py, "pos"), pos + 1)?;
        fields.set_item(intern!(py, "ref"), char::from(base))?;
        fields.set_item(intern!(py, "depth"), self.depth)?;
        fields.set_item(intern!(py, "no_calls"), ledger.no_calls[offset])?;
        fields.set_item(intern!(py, "ref_reads"), self.ref_fwd + self.ref_rev)?;
        fields.set_item(intern!(py, "ref_fwd"), self.ref_fwd)?;
        fields.set_item(intern!(py, "ref_rev"), self.ref_rev)?;
        if let Some(counts) = ledger.alleles.remove(&pos) {
            let mut alleles: Vec<(AlleleKey, Strands)> = counts.into_iter().collect();
            alleles.sort_by(|(a, [a_fwd, a_rev]), (b, [b_fwd, b_rev])| {
                (b_fwd + b_rev).cmp(&(a_fwd + a_rev)).then_with(|| a.cmp(b))
            });
            let text = |bytes: &[u8]| String::from_utf8_lossy(bytes).into_owned();
            let refs = alleles.iter().map(|((reference, _), _)| text(reference));
            let alts = alleles.iter().map(|((_, alternate), _)| text(alternate));
            let count = |pick: fn(&Strands) -> i64| -> PyResult<Bound<'_, PyTuple>> {
                PyTuple::new(py, alleles.iter().map(|(_, counts)| pick(counts)))
            };
            let refs = PyTuple::new(py, refs.collect::<Vec<_>>())?;
            fields.set_item(intern!(py, "alt_refs"), refs)?;
            let alts = PyTuple::new(py, alts.collect::<Vec<_>>())?;
            fields.set_item(intern!(py, "alts"), alts)?;
            fields.set_item(intern!(py, "alt_reads"), count(|[fwd, rev]| fwd + rev)?)?;
            fields.set_item(intern!(py, "alt_fwd"), count(|[fwd, _]| *fwd)?)?;
            fields.set_item(intern!(py, "alt_rev"), count(|[_, rev]| *rev)?)?;
        }
        let class = TABULATED_BASE.get_or_try_init(py, || {
            py.import("streampile._table")?
                .getattr("TabulatedBase")
                .map(Bound::unbind)
        })?;
        class
            .bind(py)
            .call((), Some(&fields))
            .map(|row| Some(row.unbind()))
    }
}

/// One CIGAR operator of a read: its kind, its reference start, its query start, and length.
#[derive(Clone, Copy)]
struct Segment {
    kind: SegmentKind,
    reference: i64,
    query: i64,
    length: i64,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum SegmentKind {
    Aligned,
    Inserted,
    Deleted,
    Skipped,
    Clipped,
}

/// A difference from the reference within one read, in reference and query coordinates.
#[derive(Clone, Copy)]
struct Event {
    ref_start: i64,
    ref_end: i64,
    query_start: i64,
    query_end: i64,
    is_indel: bool,
}

/// Adjacent differences of one read, and whether a matching aligned base borders each end.
///
/// `floor` is the lowest position the run's allele may be left-aligned to: the end of the read's
/// previous difference or reference skip, or else the read's first aligned base.
struct Run {
    events: Vec<Event>,
    anchored: bool,
    closed: bool,
    floor: i64,
}

/// Collects the runs of one read as its aligned bases and indels are walked in order.
struct Runs {
    finished: Vec<Run>,
    events: Vec<Event>,
    anchored: bool,
    floor: i64,
}

impl Runs {
    fn border(&mut self, aligned: bool) {
        if let Some(last) = self.events.last() {
            let floor = last.ref_end;
            self.finished.push(Run {
                events: std::mem::take(&mut self.events),
                anchored: self.anchored,
                closed: aligned,
                floor: self.floor,
            });
            self.floor = floor;
        }
        self.anchored = aligned;
    }
}

/// A normalized VCF allele seen in one read, and the reference stretch it accounts for.
struct Allele {
    pos: i64,
    reference: Vec<u8>,
    alternate: Vec<u8>,
    start: i64,
    end: i64,
}

/// What one read is counted for: its alleles, and the stretches of the alleles it is not
/// counted for and of those that spell the reference.
#[derive(Default)]
struct Counted {
    alleles: Vec<Allele>,
    dropped: Vec<(i64, i64)>,
    matched: Vec<(i64, i64)>,
}

/// One read, decoded as `tabulate` reads it.
#[derive(Default)]
struct Read {
    start: i64,
    reverse: bool,
    segments: Vec<Segment>,
    sequence: Vec<u8>,
    qualities: Vec<u8>,
}

impl Read {
    fn fill(&mut self, record: &bam::Record) -> PyResult<()> {
        let invalid = |error: std::io::Error| PyValueError::new_err(error.to_string());
        self.start = record
            .alignment_start()
            .transpose()
            .map_err(invalid)?
            .map_or(-1, |start| usize::from(start) as i64 - 1);
        self.reverse = record.flags().is_reverse_complemented();
        self.sequence.clear();
        self.sequence.extend(record.sequence().iter());
        self.qualities.clear();
        self.qualities
            .extend_from_slice(record.quality_scores().as_bytes());
        self.segments.clear();
        let (mut reference, mut query) = (self.start, 0);
        for op in record.cigar().iter() {
            let op = op.map_err(invalid)?;
            let length = op.len() as i64;
            let segment = |kind| Segment {
                kind,
                reference,
                query,
                length,
            };
            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    self.segments.push(segment(SegmentKind::Aligned));
                    reference += length;
                    query += length;
                }
                Kind::Insertion => {
                    self.segments.push(segment(SegmentKind::Inserted));
                    query += length;
                }
                Kind::Deletion => {
                    self.segments.push(segment(SegmentKind::Deleted));
                    reference += length;
                }
                Kind::Skip => {
                    self.segments.push(segment(SegmentKind::Skipped));
                    reference += length;
                }
                Kind::SoftClip => {
                    self.segments.push(segment(SegmentKind::Clipped));
                    query += length;
                }
                Kind::HardClip | Kind::Pad => {}
            }
        }
        Ok(())
    }

    /// The runs of adjacent differences from the reference, in alignment order.
    fn runs(&self, py: Python<'_>, reference: &mut Reference) -> PyResult<Vec<Run>> {
        let mut runs = Runs {
            finished: Vec::new(),
            events: Vec::new(),
            anchored: false,
            floor: self.start,
        };
        for segment in &self.segments {
            match segment.kind {
                SegmentKind::Aligned => {
                    let bases =
                        reference.get(py, segment.reference, segment.reference + segment.length)?;
                    let read = slice(
                        &self.sequence,
                        segment.query,
                        segment.query + segment.length,
                    );
                    let mut previous = -1;
                    if read != bases {
                        for (offset, (base, ref_base)) in read.iter().zip(bases).enumerate() {
                            if base != ref_base {
                                let offset = offset as i64;
                                if offset > previous + 1 {
                                    runs.border(true);
                                }
                                let at = segment.reference + offset;
                                let query_at = segment.query + offset;
                                runs.events.push(Event {
                                    ref_start: at,
                                    ref_end: at + 1,
                                    query_start: query_at,
                                    query_end: query_at + 1,
                                    is_indel: false,
                                });
                                previous = offset;
                            }
                        }
                    }
                    let aligned = bases.len() as i64;
                    if previous < aligned - 1 {
                        runs.border(true);
                    }
                    if aligned < segment.length {
                        runs.border(false);
                    }
                }
                SegmentKind::Inserted => runs.events.push(Event {
                    ref_start: segment.reference,
                    ref_end: segment.reference,
                    query_start: segment.query,
                    query_end: segment.query + segment.length,
                    is_indel: true,
                }),
                SegmentKind::Deleted => runs.events.push(Event {
                    ref_start: segment.reference,
                    ref_end: segment.reference + segment.length,
                    query_start: segment.query,
                    query_end: segment.query,
                    is_indel: true,
                }),
                SegmentKind::Skipped | SegmentKind::Clipped => {
                    runs.border(false);
                    if segment.kind == SegmentKind::Skipped {
                        runs.floor = segment.reference + segment.length;
                    }
                }
            }
        }
        runs.border(false);
        Ok(runs.finished)
    }
}

/// The options of a tabulation, and the reference windows of each contig it has read.
struct Options {
    min_base_quality: u8,
    min_mapping_quality: u8,
    exclude_flags: u16,
}

impl Options {
    fn accepts(&self, record: &bam::Record) -> bool {
        let flags = record.flags();
        !flags.is_unmapped()
            && flags.bits() & self.exclude_flags == 0
            && record.mapping_quality().map_or(255, u8::from) >= self.min_mapping_quality
            && !record.sequence().is_empty()
            && !record.cigar().is_empty()
    }

    /// The alleles a read is counted for, and the stretches of the alleles it is not counted
    /// for and of those that spell the reference.
    fn alleles(&self, py: Python<'_>, read: &Read, reference: &mut Reference) -> PyResult<Counted> {
        let qualities =
            (self.min_base_quality > 0 && !read.qualities.is_empty()).then_some(&read.qualities);
        let mut counted = Counted::default();
        for run in read.runs(py, reference)? {
            let (first, last) = (run.events[0], run.events[run.events.len() - 1]);
            let (mut ref_start, mut query_start) = (first.ref_start, first.query_start);
            let (ref_end, query_end) = (last.ref_end, last.query_end);
            if first.is_indel {
                if !run.anchored {
                    counted.dropped.push((ref_start, ref_end));
                    continue;
                }
                ref_start -= 1;
                query_start -= 1;
            }
            let reference_bases = reference.get(py, ref_start, ref_end)?.to_vec();
            let alternate = python_slice(&read.sequence, query_start, query_end).to_vec();
            let opens_with_insertion = first.is_indel && first.ref_end == first.ref_start;
            let closes_with_deletion = last.is_indel && last.ref_end > last.ref_start;
            let judged_from = if opens_with_insertion {
                query_start
            } else {
                first.query_start
            };
            let judged_to = if closes_with_deletion {
                query_end + 1
            } else {
                query_end
            };
            let under_floor = || -> PyResult<bool> {
                let Some(qualities) = qualities else {
                    return Ok(false);
                };
                let judged = python_slice(qualities, judged_from, judged_to);
                let lowest = judged
                    .iter()
                    .min()
                    .ok_or_else(|| PyValueError::new_err("min() iterable argument is empty"))?;
                Ok(*lowest < self.min_base_quality)
            };
            if (closes_with_deletion && !run.closed)
                || under_floor()?
                || alternate.contains(&b'N')
                || !reference_bases.iter().copied().all(is_acgt)
            {
                counted.dropped.push((ref_start, ref_end));
                continue;
            }
            if reference_bases == alternate {
                counted.matched.push((ref_start, ref_end));
                continue;
            }
            let normalized = normalize(
                ref_start,
                reference_bases,
                alternate,
                |start, end| reference.get(py, start, end).map(<[u8]>::to_vec),
                run.floor,
            )?;
            let Some((pos, reference_bases, alternate)) = normalized else {
                counted.dropped.push((run.floor, ref_end));
                continue;
            };
            let start = ref_start.min(pos);
            let end = ref_end.max(pos + reference_bases.len() as i64);
            if !reference
                .get(py, start, ref_start)?
                .iter()
                .copied()
                .all(is_acgt)
            {
                counted.dropped.push((start, end));
                continue;
            }
            counted.alleles.push(Allele {
                pos,
                reference: reference_bases,
                alternate,
                start,
                end,
            });
        }
        Ok(counted)
    }

    /// The aligned positions of a read under the floor or over a non-ACGT base, and its `N`s.
    fn uninformative(
        &self,
        py: Python<'_>,
        read: &Read,
        reference: &mut Reference,
    ) -> PyResult<(Vec<i64>, Vec<i64>)> {
        let floor = self.min_base_quality;
        let low = floor > 0 && !read.qualities.is_empty();
        let mut uninformative = Vec::new();
        let mut no_calls = Vec::new();
        for segment in &read.segments {
            if segment.kind != SegmentKind::Aligned {
                continue;
            }
            let end = segment.query + segment.length;
            let at = |offset: i64| segment.reference + offset;
            let from = segment.query.clamp(0, read.sequence.len() as i64);
            for (offset, &base) in python_slice(&read.sequence, segment.query, end)
                .iter()
                .enumerate()
            {
                if base == b'N' {
                    no_calls.push(at(from + offset as i64 - segment.query));
                }
            }
            if low {
                let from = segment.query.clamp(0, read.qualities.len() as i64);
                for (offset, &quality) in python_slice(&read.qualities, segment.query, end)
                    .iter()
                    .enumerate()
                {
                    if quality < floor {
                        uninformative.push(at(from + offset as i64 - segment.query));
                    }
                }
            }
            let bases = reference.get(py, segment.reference, segment.reference + segment.length)?;
            for (offset, &base) in bases.iter().enumerate() {
                if !is_acgt(base) {
                    uninformative.push(segment.reference + offset as i64);
                }
            }
        }
        Ok((uninformative, no_calls))
    }
}

/// A Python slice of bytes from `start` to `end`, with negative bounds counted from the end.
fn python_slice(bytes: &[u8], start: i64, end: i64) -> &[u8] {
    let length = bytes.len() as i64;
    let bound = |value: i64| if value < 0 { value + length } else { value };
    slice(bytes, bound(start), bound(end))
}

/// Trims and left-aligns an allele as `bcftools norm` does, or returns `None` for an allele
/// that changes nothing or would move past `floor`.
fn normalize<T: Copy + Eq>(
    mut pos: i64,
    mut reference: Vec<T>,
    mut alternate: Vec<T>,
    mut bases: impl FnMut(i64, i64) -> PyResult<Vec<T>>,
    floor: i64,
) -> PyResult<Option<Normalized<T>>> {
    if reference == alternate {
        return Ok(None);
    }
    loop {
        let (Some(last), Some(other)) = (reference.last(), alternate.last()) else {
            return Err(PyIndexError::new_err("string index out of range"));
        };
        if last != other {
            break;
        }
        reference.pop();
        alternate.pop();
        if reference.is_empty() || alternate.is_empty() {
            if pos - 1 < floor {
                return Ok(None);
            }
            pos -= 1;
            let base = bases(pos, pos + 1)?;
            reference.splice(0..0, base.iter().copied());
            alternate.splice(0..0, base);
        }
    }
    let trimmed = reference
        .iter()
        .zip(&alternate)
        .take(reference.len().min(alternate.len()).saturating_sub(1))
        .take_while(|(a, b)| a == b)
        .count();
    reference.drain(..trimmed);
    alternate.drain(..trimmed);
    Ok(Some((pos + trimmed as i64, reference, alternate)))
}

/// Trim and left-align an allele as `bcftools norm` does.
#[pyfunction]
#[pyo3(signature = (pos, reference_bases, alternate_bases, reference, floor = 0))]
pub(crate) fn normalize_allele(
    pos: i64,
    reference_bases: &str,
    alternate_bases: &str,
    reference: &Bound<'_, PyAny>,
    floor: i64,
) -> PyResult<Option<(i64, String, String)>> {
    let normalized = normalize(
        pos,
        reference_bases.chars().collect(),
        alternate_bases.chars().collect(),
        |start, end| {
            let bases = reference.call1((start, end))?;
            Ok(bases.cast::<PyString>()?.to_str()?.chars().collect())
        },
        floor,
    )?;
    Ok(normalized.map(|(pos, reference, alternate)| {
        (
            pos,
            reference.into_iter().collect(),
            alternate.into_iter().collect(),
        )
    }))
}

/// The per-read and per-base work of a `Tabulator`, over one reference.
#[pyclass(module = "streampile._native")]
pub(crate) struct Tabulation {
    fasta: Py<PyAny>,
    options: Options,
    contigs: HashMap<String, Reference>,
}

#[pymethods]
impl Tabulation {
    #[new]
    #[pyo3(signature = (
        reference,
        *,
        min_base_quality = Int::from(0),
        min_mapping_quality = Int::from(0),
        exclude_flags = Int::from(i64::from(DEFAULT_EXCLUDE_FLAGS.bits())),
    ))]
    #[pyo3(
        text_signature = "(reference, *, min_base_quality=0, min_mapping_quality=0, \
        exclude_flags=3840)"
    )]
    fn new(
        reference: &Bound<'_, PyAny>,
        min_base_quality: Int,
        min_mapping_quality: Int,
        exclude_flags: Int,
    ) -> PyResult<Self> {
        Ok(Self {
            fasta: reference.clone().unbind(),
            options: Options {
                min_base_quality: min_base_quality.of("min_base_quality", u8::MAX)?,
                min_mapping_quality: min_mapping_quality.of("min_mapping_quality", u8::MAX)?,
                exclude_flags: exclude_flags.of("exclude_flags", u16::MAX)?,
            },
            contigs: HashMap::new(),
        })
    }

    /// Whether a read passes the flag and mapping-quality filters and has bases to count.
    fn accepts(&self, record: &Bound<'_, PyAny>) -> PyResult<bool> {
        let mut bam = bam::Record::default();
        bridge::read(record, &mut bam)?;
        Ok(self.options.accepts(&bam))
    }

    /// The alleles a read is counted for, as `(pos, ref, alt, start, end)`, and the stretches of
    /// those it is not counted for.
    #[allow(clippy::type_complexity)]
    fn alleles(
        &mut self,
        py: Python<'_>,
        record: &Bound<'_, PyAny>,
    ) -> PyResult<(Vec<(i64, String, String, i64, i64)>, Vec<(i64, i64)>)> {
        let contig = record.getattr("reference_name")?;
        if contig.is_none() {
            let name: Option<String> = record.getattr("query_name")?.extract()?;
            return Err(PyValueError::new_err(format!(
                "Read {} is not mapped.",
                name.as_deref().unwrap_or("None")
            )));
        }
        let contig = contig.cast_into::<PyString>()?;
        let key = contig.to_str()?.to_owned();
        let reference = match self.contigs.entry(key) {
            std::collections::hash_map::Entry::Occupied(entry) => entry.into_mut(),
            std::collections::hash_map::Entry::Vacant(entry) => {
                entry.insert(Reference::new(self.fasta.bind(py), &contig)?)
            }
        };
        let mut bam = bam::Record::default();
        bridge::read(record, &mut bam)?;
        let mut read = Read::default();
        read.fill(&bam)?;
        let counted = self.options.alleles(py, &read, reference)?;
        let text = |bytes: Vec<u8>| String::from_utf8_lossy(&bytes).into_owned();
        Ok((
            counted
                .alleles
                .into_iter()
                .map(|allele| {
                    (
                        allele.pos,
                        text(allele.reference),
                        text(allele.alternate),
                        allele.start,
                        allele.end,
                    )
                })
                .collect(),
            counted.dropped,
        ))
    }

    /// The rows of every base of the spans of one contig, in order, from reads fetched span by
    /// span from an indexed alignment file.
    fn tabulate_contig(
        &self,
        py: Python<'_>,
        alignments: &Bound<'_, PyAny>,
        contig: &Bound<'_, PyString>,
        spans: Vec<(i64, i64)>,
    ) -> PyResult<ContigRows> {
        Ok(ContigRows {
            alignments: alignments.clone().unbind(),
            contig: contig.clone().unbind(),
            reference: Reference::new(self.fasta.bind(py), contig)?,
            options: Options {
                min_base_quality: self.options.min_base_quality,
                min_mapping_quality: self.options.min_mapping_quality,
                exclude_flags: self.options.exclude_flags,
            },
            starts: spans.iter().map(|&(start, _)| start).collect(),
            spans,
            ledgers: BTreeMap::new(),
            next_span: 0,
            fetched_to: -1,
            rows: None,
            record: bam::Record::default(),
            read: Read::default(),
        })
    }
}

/// The rows of every base of the spans of one contig, made one at a time.
#[pyclass(module = "streampile._native")]
pub(crate) struct ContigRows {
    alignments: Py<PyAny>,
    contig: Py<PyString>,
    reference: Reference,
    options: Options,
    spans: Vec<(i64, i64)>,
    starts: Vec<i64>,
    ledgers: BTreeMap<usize, Ledger>,
    next_span: usize,
    fetched_to: i64,
    rows: Option<Rows>,
    record: bam::Record,
    read: Read,
}

#[pymethods]
impl ContigRows {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<Py<PyAny>>> {
        loop {
            if let Some(rows) = self.rows.as_mut()
                && let Some(row) = rows.next(py, self.contig.bind(py))?
            {
                return Ok(Some(row));
            }
            let Some(&(start, end)) = self.spans.get(self.next_span) else {
                return Ok(None);
            };
            let index = self.next_span;
            self.next_span += 1;
            let records = self
                .alignments
                .bind(py)
                .call_method1("fetch", (self.contig.bind(py), start, end))?;
            for record in records.try_iter()? {
                bridge::read(&record?, &mut self.record)?;
                self.count()?;
            }
            self.fetched_to = end;
            let ledger = self
                .ledgers
                .remove(&index)
                .unwrap_or_else(|| Ledger::new(start, end));
            self.rows = Some(Rows::new(py, ledger, &mut self.reference)?);
        }
    }
}

impl ContigRows {
    fn count(&mut self) -> PyResult<()> {
        let record = &self.record;
        let start = record
            .alignment_start()
            .transpose()
            .map_err(|error| PyValueError::new_err(error.to_string()))?
            .map_or(-1, |start| usize::from(start) as i64 - 1);
        if start < self.fetched_to || !self.options.accepts(record) {
            return Ok(());
        }
        let end = end_of(record, start);
        let touched = self.touch(start, end);
        if touched.is_empty() {
            return Ok(());
        }
        Python::attach(|py| {
            self.read.fill(&self.record)?;
            let read = &self.read;
            let reverse = read.reverse;
            let counted = self.options.alleles(py, read, &mut self.reference)?;
            let (uninformative, no_calls) =
                self.options.uninformative(py, read, &mut self.reference)?;
            let mut excluded = uninformative;
            excluded.extend_from_slice(&no_calls);
            let ledgers = &mut self.ledgers;
            let mut each = |apply: &mut dyn FnMut(&mut Ledger)| {
                for index in &touched {
                    if let Some(ledger) = ledgers.get_mut(index) {
                        apply(ledger);
                    }
                }
            };
            for &pos in &no_calls {
                each(&mut |ledger| ledger.add_no_call(pos));
            }
            for &(start, end) in &counted.dropped {
                excluded.extend(start..end);
            }
            for allele in &counted.alleles {
                let after = allele.pos + allele.reference.len() as i64;
                each(&mut |ledger| {
                    ledger.add(allele.start, allele.pos, true, reverse);
                    ledger.add(allele.pos, after, false, reverse);
                    ledger.add(after, allele.end, true, reverse);
                    ledger.add_allele(allele, reverse);
                });
                excluded.extend(allele.start..allele.end);
            }
            for &(start, end) in &counted.matched {
                each(&mut |ledger| ledger.add(start, end, true, reverse));
                excluded.extend(start..end);
            }
            excluded.sort_unstable();
            excluded.dedup();
            for segment in &read.segments {
                if segment.kind != SegmentKind::Aligned {
                    continue;
                }
                let (block_start, block_end) =
                    (segment.reference, segment.reference + segment.length);
                let from = excluded.partition_point(|&pos| pos < block_start);
                let to = excluded.partition_point(|&pos| pos < block_end);
                let mut cursor = block_start;
                for &pos in &excluded[from..to] {
                    each(&mut |ledger| ledger.add(cursor, pos, true, reverse));
                    cursor = pos + 1;
                }
                each(&mut |ledger| ledger.add(cursor, block_end, true, reverse));
            }
            Ok(())
        })
    }

    /// The indices of the spans a read from `start` to `end` overlaps, with their ledgers opened
    /// as needed.
    fn touch(&mut self, start: i64, end: i64) -> Vec<usize> {
        let mut touched = Vec::new();
        let mut index = self
            .starts
            .partition_point(|&first| first <= start)
            .saturating_sub(1);
        while let Some(&(span_start, span_end)) = self.spans.get(index) {
            if span_start >= end {
                break;
            }
            if span_end > start {
                self.ledgers
                    .entry(index)
                    .or_insert_with(|| Ledger::new(span_start, span_end));
                touched.push(index);
            }
            index += 1;
        }
        touched
    }
}

/// Where a read ends on the reference, one past its start for one that spans no base, as pysam's
/// `reference_end` gives.
fn end_of(record: &bam::Record, start: i64) -> i64 {
    use noodles::sam::alignment::record::Cigar as _;
    let cigar = record.cigar();
    if record.flags().is_unmapped() || cigar.is_empty() {
        start
    } else {
        start + (cigar.alignment_span().unwrap_or(0) as i64).max(1)
    }
}
