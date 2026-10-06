//! `Pileup` and `PileupRead`, owned snapshots of a column that outlive the builder's next move.

use std::sync::Arc;

use noodles::bam;
use noodles::sam::alignment::record::Cigar as _;
use noodles::sam::alignment::record::cigar::op::Kind;
use pyo3::exceptions::{PyIndexError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::pyclass::CompareOp;
use pyo3::sync::PyOnceLock;
use pyo3::types::{PyDict, PyList, PyString, PyTuple, PyType};

use super::builder::Bridged;
use super::{bridge, to_python};
use crate::builder::{MateCigar, other_end};
use crate::footprint::{Footprint, Located};
use crate::pileup::{
    EntryKind, MISSING_BASE_QUALITY, NONE, OtherEnd, Template, name_hash, name_of,
    template_end_distance, templates_kept,
};
use crate::{DEFAULT_MIN_BASE_QUALITY, PileupEntry};

const ABSENT: i64 = i64::MIN;

const FIELDS: [&str; 6] = [
    "alignment",
    "query_position",
    "query_position_or_next",
    "pileup_type",
    "insertion_offset",
    "insertion_length",
];

static PILEUP_READ_TYPES: PyOnceLock<[Py<PyAny>; 4]> = PyOnceLock::new();

/// One read at one position: the record, the fields of a `PileupRead`, and what is worked out
/// about the record once, where a builder made it.
#[derive(Clone)]
pub(crate) struct Held {
    record: Arc<Bridged>,
    kind: EntryKind,
    query_position: i64,
    query_position_or_next: i64,
    insertion_offset: i64,
    insertion_length: i64,
    other_end: Option<OtherEnd>,
    query_length: Option<u32>,
    name_hash: u64,
}

impl Held {
    /// The entry of a builder's pileup.
    pub(crate) fn of(entry: &PileupEntry<'_, super::builder::PyRecord>) -> Self {
        let present = |value: Option<usize>| value.map_or(ABSENT, |value| value as i64);
        let live = entry.live;
        Self {
            record: Arc::clone(entry.source_record().shared()),
            kind: entry.kind(),
            query_position: present(entry.query_position()),
            query_position_or_next: present(entry.query_position_or_next()),
            insertion_offset: present(entry.insertion_offset()),
            insertion_length: entry.insertion_length() as i64,
            other_end: Some(live.other_end),
            query_length: Some(live.footprint.query_length),
            name_hash: live.name_hash,
        }
    }

    /// An entry of a read placed by its footprint, as a builder would make it at a position.
    fn located(record: &Arc<Bridged>, kind: EntryKind, offset: u32, length: u32) -> Self {
        let offset = if offset == NONE {
            ABSENT
        } else {
            i64::from(offset)
        };
        Self {
            record: Arc::clone(record),
            kind,
            query_position: if kind == EntryKind::Base {
                offset
            } else {
                ABSENT
            },
            query_position_or_next: if matches!(kind, EntryKind::Base | EntryKind::Deletion) {
                offset
            } else {
                ABSENT
            },
            insertion_offset: if kind == EntryKind::Insertion {
                offset
            } else {
                ABSENT
            },
            insertion_length: if kind == EntryKind::Insertion {
                i64::from(length)
            } else {
                0
            },
            other_end: None,
            query_length: None,
            name_hash: name_hash(&record.record),
        }
    }

    fn bam(&self) -> &bam::Record {
        &self.record.record
    }

    fn base(&self) -> PyResult<Option<u8>> {
        let sequence = self.bam().sequence();
        if self.query_position == ABSENT || sequence.is_empty() {
            return Ok(None);
        }
        let index = index(self.query_position, sequence.len())
            .ok_or_else(|| PyIndexError::new_err("string index out of range"))?;
        Ok(sequence.get(index))
    }

    fn quality(&self) -> PyResult<Option<u8>> {
        let offset = self.query_position_or_next;
        if offset == ABSENT {
            return Ok(None);
        }
        let record = self.bam();
        let qualities = record.quality_scores().as_bytes();
        if qualities.is_empty() {
            return Ok((offset < record.sequence().len() as i64).then_some(MISSING_BASE_QUALITY));
        }
        let index = index(offset, qualities.len())
            .ok_or_else(|| PyIndexError::new_err("array index out of range"))?;
        Ok(Some(qualities[index]))
    }

    fn passes(&self, floor: i64) -> PyResult<bool> {
        Ok(self
            .quality()?
            .is_some_and(|quality| i64::from(quality) >= floor))
    }

    fn inserted_bases(&self) -> Option<String> {
        let sequence = self.bam().sequence();
        if self.insertion_offset == ABSENT || sequence.is_empty() {
            return None;
        }
        let (start, end) = slice(
            self.insertion_offset,
            self.insertion_offset.saturating_add(self.insertion_length),
            sequence.len(),
        );
        Some(
            (start..end)
                .filter_map(|index| sequence.get(index).map(char::from))
                .collect(),
        )
    }

    fn inserted_qualities(&self) -> Option<Vec<u8>> {
        if self.insertion_offset == ABSENT {
            return None;
        }
        let record = self.bam();
        let qualities = record.quality_scores().as_bytes();
        if qualities.is_empty() {
            let stored = record.sequence().len();
            let length = usize::try_from(self.insertion_length).unwrap_or(0);
            return (stored > 0).then(|| vec![MISSING_BASE_QUALITY; length]);
        }
        let (start, end) = slice(
            self.insertion_offset,
            self.insertion_offset.saturating_add(self.insertion_length),
            qualities.len(),
        );
        Some(qualities[start..end].to_vec())
    }

    fn five_prime_distance(&self) -> Option<i64> {
        let offset = self.query_position;
        if offset < 0 {
            return None;
        }
        let record = self.bam();
        if !record.flags().is_reverse_complemented() {
            return Some(offset);
        }
        let length = match self.query_length {
            Some(length) => i64::from(length),
            None => i64::try_from(record.cigar().read_length().ok()?).ok()?,
        };
        Some(length - offset - 1).filter(|&distance| distance >= 0)
    }

    fn template_end_distance(&self, position: Option<i64>) -> PyResult<Option<i64>> {
        let record = self.bam();
        let Some(position) = position.or_else(|| self.position_of_base()) else {
            return Ok(None);
        };
        let other_end = match self.other_end {
            Some(other_end) => other_end,
            None => other_end_of(record).map_err(to_python)?,
        };
        let distance = template_end_distance(record, other_end, position).map_err(to_python)?;
        Ok(distance.map(|distance| distance as i64))
    }

    /// The reference position of the read's base at the query position, for an entry made by
    /// hand, which knows no position of its own.
    fn position_of_base(&self) -> Option<i64> {
        if self.kind != EntryKind::Base || self.query_position < 0 {
            return None;
        }
        let record = self.bam();
        let mut reference = usize::from(record.alignment_start()?.ok()?) as i64 - 1;
        let mut query = 0;
        for op in record.cigar().iter() {
            let op = op.ok()?;
            let length = op.len() as i64;
            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    if (query..query + length).contains(&self.query_position) {
                        return Some(reference + self.query_position - query);
                    }
                    reference += length;
                    query += length;
                }
                Kind::Insertion | Kind::SoftClip => query += length,
                Kind::Deletion | Kind::Skip => reference += length,
                Kind::HardClip | Kind::Pad => {}
            }
        }
        None
    }

    fn template(&self) -> Template<'_> {
        Template {
            hash: self.name_hash,
            name: name_of(self.bam()),
        }
    }

    fn fields<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let optional = |value: i64| (value != ABSENT).then_some(value);
        (
            self.record.object.clone_ref(py),
            optional(self.query_position),
            optional(self.query_position_or_next),
            pileup_read_type(py, self.kind)?,
            optional(self.insertion_offset),
            self.insertion_length,
        )
            .into_pyobject(py)
    }
}

/// Where the template of a read made by hand ends, worked out from the read alone.
fn other_end_of(record: &bam::Record) -> crate::Result<OtherEnd> {
    let (Some(reference_id), Some(start)) = (
        record.reference_sequence_id().transpose()?,
        record.alignment_start().transpose()?,
    ) else {
        return Ok(OtherEnd::None);
    };
    let start = usize::from(start) as i64 - 1;
    let end = start + record.cigar().alignment_span()? as i64;
    other_end(record, reference_id, start, end, MateCigar::Unsearched)
}

/// The index a Python index of a sequence of this length refers to, if any.
fn index(position: i64, length: usize) -> Option<usize> {
    let length = length as i64;
    let position = if position < 0 {
        position + length
    } else {
        position
    };
    (0..length).contains(&position).then_some(position as usize)
}

/// The bounds a Python slice from `start` to `stop` of a sequence of this length refers to.
fn slice(start: i64, stop: i64, length: usize) -> (usize, usize) {
    let bound = |value: i64| {
        let value = if value < 0 {
            value + length as i64
        } else {
            value
        };
        value.clamp(0, length as i64) as usize
    };
    let (start, stop) = (bound(start), bound(stop));
    (start, stop.max(start))
}

fn kind_of(name: &str) -> PyResult<EntryKind> {
    Ok(match name {
        "base" => EntryKind::Base,
        "deletion" => EntryKind::Deletion,
        "insertion" => EntryKind::Insertion,
        "skip" => EntryKind::Skip,
        _ => {
            return Err(PyValueError::new_err(format!(
                "'{name}' is not a valid PileupReadType"
            )));
        }
    })
}

/// The `PileupReadType` member of an entry kind.
fn pileup_read_type(py: Python<'_>, kind: EntryKind) -> PyResult<Bound<'_, PyAny>> {
    let members = PILEUP_READ_TYPES.get_or_try_init(py, || -> PyResult<_> {
        let types = py.import("streampile._pileup")?.getattr("PileupReadType")?;
        let member = |name: &str| types.getattr(name).map(Bound::unbind);
        Ok([
            member("base")?,
            member("deletion")?,
            member("insertion")?,
            member("skip")?,
        ])
    })?;
    let index = match kind {
        EntryKind::Base => 0,
        EntryKind::Deletion => 1,
        EntryKind::Insertion => 2,
        EntryKind::Skip => 3,
    };
    Ok(members[index].bind(py).clone())
}

fn optional(value: Option<i64>) -> i64 {
    value.unwrap_or(ABSENT)
}

/// One read at one pileup position.
///
/// A read holding a base, a deletion, or a reference skip (the CIGAR `N` operator) at a position
/// appears once. A read with an insertion right after the position appears again as an insertion
/// entry, as does a read whose alignment opens with an insertion, at the position before its first
/// aligned base. So an insertion at either end of an alignment is reported: htslib reports only one
/// that closes an alignment, and fgbio only one that opens it, at offset 0 even after a soft clip.
///
/// A skip entry holds no base, no quality, and no query offset. htslib flags the same entry as
/// both `is_del` and `is_refskip` and gives it the offset and quality of the read's next base;
/// here a skip is not a deletion, and has no quality to pass a floor with.
///
/// A `PileupRead` unpacks, indexes, compares, and hashes as the tuple of its six fields, as a
/// `NamedTuple` does. Its bases and qualities are those of the read as it was when it was piled
/// up, while `alignment` is the very object given to the builder, so it can be changed, e.g.
/// tagged, and written on.
///
/// Attributes:
///     alignment: the read.
///     query_position: the 0-based query offset of the read's base at the position, or `None`
///         for a deletion, a skip, or an insertion.
///     query_position_or_next: the query offset of the read's base at the position, or of the
///         read's next base for a deletion, as htslib reports it, or `None` for a skip, an
///         insertion, or a deletion no base follows.
///     pileup_type: whether the read holds a base, a deletion, a skip, or an insertion.
///     insertion_offset: the query offset of the first inserted base of an insertion entry.
///     insertion_length: the number of inserted bases of an insertion entry.
#[pyclass(module = "streampile", name = "PileupRead", frozen)]
pub(crate) struct PileupRead {
    held: Held,
    position: Option<i64>,
}

#[pymethods]
impl PileupRead {
    #[new]
    #[pyo3(signature = (
        alignment,
        query_position,
        query_position_or_next,
        pileup_type,
        insertion_offset = None,
        insertion_length = 0,
    ))]
    fn new(
        alignment: &Bound<'_, PyAny>,
        query_position: Option<i64>,
        query_position_or_next: Option<i64>,
        pileup_type: &str,
        insertion_offset: Option<i64>,
        insertion_length: i64,
    ) -> PyResult<Self> {
        let mut record = bam::Record::default();
        bridge::read(alignment, &mut record)?;
        let name_hash = name_hash(&record);
        Ok(Self {
            held: Held {
                record: Arc::new(Bridged {
                    object: alignment.clone().unbind(),
                    record,
                }),
                kind: kind_of(pileup_type)?,
                query_position: optional(query_position),
                query_position_or_next: optional(query_position_or_next),
                insertion_offset: optional(insertion_offset),
                insertion_length,
                other_end: None,
                query_length: None,
                name_hash,
            },
            position: None,
        })
    }

    #[classattr]
    #[pyo3(name = "_fields")]
    fn field_names(py: Python<'_>) -> PyResult<Bound<'_, PyTuple>> {
        PyTuple::new(py, FIELDS)
    }

    /// The read.
    #[getter]
    fn alignment(&self, py: Python<'_>) -> Py<PyAny> {
        self.held.record.object.clone_ref(py)
    }

    /// The 0-based query offset of the read's base at the position, or `None` without one.
    #[getter]
    fn query_position(&self) -> Option<i64> {
        (self.held.query_position != ABSENT).then_some(self.held.query_position)
    }

    /// The query offset of the read's base here, or of its next base for a deletion, or `None`.
    #[getter]
    fn query_position_or_next(&self) -> Option<i64> {
        (self.held.query_position_or_next != ABSENT).then_some(self.held.query_position_or_next)
    }

    /// Whether the read holds a base, a deletion, a skip, or an insertion.
    #[getter]
    fn pileup_type<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        pileup_read_type(py, self.held.kind)
    }

    /// The query offset of the first inserted base of an insertion entry.
    #[getter]
    fn insertion_offset(&self) -> Option<i64> {
        (self.held.insertion_offset != ABSENT).then_some(self.held.insertion_offset)
    }

    /// The number of inserted bases of an insertion entry.
    #[getter]
    fn insertion_length(&self) -> i64 {
        self.held.insertion_length
    }

    /// The upper-cased read base at the position, or `None` without one.
    #[getter]
    fn base(&self) -> PyResult<Option<char>> {
        Ok(self.held.base()?.map(char::from))
    }

    /// The base quality at the position, or of the read's next base for a deletion.
    ///
    /// A read with bases but no stored qualities (QUAL `*`) has quality 255 at every base, as in
    /// htslib, so it passes every floor. It is `None` where there is no base to take it from: for
    /// an insertion, a deletion no base follows, or a read with no stored bases (SEQ `*`).
    #[getter]
    fn qual(&self) -> PyResult<Option<u8>> {
        self.held.quality()
    }

    /// Whether the read has a deletion at the position.
    #[getter]
    fn is_del(&self) -> bool {
        self.held.kind == EntryKind::Deletion
    }

    /// Whether the read holds an `N` base, a no-call, at the position.
    #[getter]
    fn is_no_call(&self) -> PyResult<bool> {
        Ok(self.held.base()? == Some(b'N'))
    }

    /// Whether this entry is an insertion right after the position.
    #[getter]
    fn is_ins(&self) -> bool {
        self.held.kind == EntryKind::Insertion
    }

    /// Whether the read skips over the position, with the CIGAR `N` operator.
    #[getter]
    fn is_refskip(&self) -> bool {
        self.held.kind == EntryKind::Skip
    }

    /// The upper-cased inserted bases of an insertion entry, or `None` for any other.
    #[getter]
    fn inserted_bases(&self) -> Option<String> {
        self.held.inserted_bases()
    }

    /// The base qualities of the inserted bases of an insertion entry, or `None`.
    ///
    /// They are 255 for a read with no stored qualities, and `None` for one with no stored bases.
    #[getter]
    fn inserted_qualities<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyList>>> {
        self.held
            .inserted_qualities()
            .map(|qualities| PyList::new(py, qualities))
            .transpose()
    }

    /// The distance of the read's base from its 5′ end, in bases as sequenced.
    ///
    /// It is the query offset for a forward read, counted from the other end for a reverse read,
    /// so soft-clipped bases count, and 0 is the first base sequenced: fgbio's
    /// `positionInReadInReadOrder` minus one. It is `None` for an entry with no base.
    #[getter]
    fn five_prime_distance(&self) -> Option<i64> {
        self.held.five_prime_distance()
    }

    /// The distance on the reference from the position to the template's other end, the
    /// unclipped 5′ end of the mate of a read in an FR pair: 0 at the mate's 5′ end.
    ///
    /// The mate's 5′ end comes from its start and its `MC` tag, counting its soft and hard clips:
    /// its start less the clips before it for a reverse read, and the end of its alignment plus
    /// the clips after it for a forward read; the template length (TLEN) is never read. It is
    /// `None` for a fragment, a read whose mate is unmapped or on another contig, a pair that is
    /// not FR, a position past the mate's 5′ end, where a read runs through its mate, and an entry
    /// made by hand that holds no base, whose position is unknown.
    ///
    /// Raises:
    ///     ValueError: for a read of an FR pair with no `MC` tag, or one that is not a CIGAR
    ///         string spanning at least one base.
    #[getter]
    fn template_end_distance(&self) -> PyResult<Option<i64>> {
        self.held.template_end_distance(self.position)
    }

    #[pyo3(name = "_asdict")]
    fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let fields = self.held.fields(py)?;
        let dict = PyDict::new(py);
        for (name, value) in FIELDS.iter().zip(fields.iter()) {
            dict.set_item(name, value)?;
        }
        Ok(dict)
    }

    #[pyo3(signature = (**changes))]
    #[pyo3(name = "_replace")]
    fn replace(&self, py: Python<'_>, changes: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
        let fields = self.asdict(py)?;
        if let Some(changes) = changes {
            for (name, value) in changes {
                if !FIELDS.contains(&name.extract::<&str>()?) {
                    return Err(PyValueError::new_err(format!(
                        "Got unexpected field names: [{}]",
                        name.repr()?
                    )));
                }
                fields.set_item(name, value)?;
            }
        }
        let read = PileupRead::new(
            &fields.as_any().get_item("alignment")?,
            fields.as_any().get_item("query_position")?.extract()?,
            fields
                .as_any()
                .get_item("query_position_or_next")?
                .extract()?,
            fields.as_any().get_item("pileup_type")?.extract()?,
            fields.as_any().get_item("insertion_offset")?.extract()?,
            fields.as_any().get_item("insertion_length")?.extract()?,
        )?;
        let unchanged =
            changes.is_none_or(|changes| !changes.contains("alignment").unwrap_or(true));
        Ok(if unchanged {
            Self {
                held: Held {
                    record: Arc::clone(&self.held.record),
                    other_end: self.held.other_end,
                    query_length: self.held.query_length,
                    ..read.held
                },
                position: self.position,
            }
        } else {
            read
        })
    }

    #[allow(clippy::unused_self)]
    fn __len__(&self) -> usize {
        FIELDS.len()
    }

    fn __getitem__<'py>(
        &self,
        py: Python<'py>,
        index: &Bound<'py, PyAny>,
    ) -> PyResult<Bound<'py, PyAny>> {
        self.held.fields(py)?.as_any().get_item(index)
    }

    fn __iter__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        Ok(self.held.fields(py)?.as_any().try_iter()?.into_any())
    }

    fn __hash__(&self, py: Python<'_>) -> PyResult<isize> {
        self.held.fields(py)?.hash()
    }

    fn __richcmp__<'py>(
        &self,
        py: Python<'py>,
        other: &Bound<'py, PyAny>,
        op: CompareOp,
    ) -> PyResult<Py<PyAny>> {
        let other = if let Ok(read) = other.cast::<PileupRead>() {
            read.get().held.fields(py)?.into_any()
        } else if other.is_instance_of::<PyTuple>() {
            other.clone()
        } else {
            return Ok(py.NotImplemented());
        };
        Ok(self.held.fields(py)?.rich_compare(other, op)?.unbind())
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let fields = self.held.fields(py)?;
        let mut parts = Vec::with_capacity(FIELDS.len());
        for (name, value) in FIELDS.iter().zip(fields.iter()) {
            parts.push(format!("{name}={}", value.repr()?));
        }
        Ok(format!("PileupRead({})", parts.join(", ")))
    }
}

/// The reads at one reference position.
///
/// Reads with an insertion right after the position are included as insertion entries, so one read
/// can have two entries. Reads that skip over the position with the CIGAR `N` operator are included
/// as skip entries. A pileup is a snapshot: it outlives the builder's next move, and its views are
/// worked out from it in Rust. It compares and hashes by its four fields.
///
/// Attributes:
///     reference_name: the name of the contig.
///     reference_pos: the 0-based position on the contig.
///     pileups: the entries of the reads at this position.
///     min_base_quality: the base quality below which bases are left out of `filtered_depth`,
///         `bases`, and `qualities`.
#[pyclass(module = "streampile", name = "Pileup", frozen)]
pub(crate) struct Pileup {
    reference_name: Py<PyString>,
    reference_pos: i64,
    min_base_quality: i64,
    entries: Vec<Held>,
    pileups: PyOnceLock<Py<PyTuple>>,
}

impl Pileup {
    /// A pileup of entries, which makes its `PileupRead`s only when asked for them.
    pub(crate) fn of(
        reference_name: Py<PyString>,
        reference_pos: i64,
        min_base_quality: i64,
        entries: Vec<Held>,
    ) -> Self {
        Self {
            reference_name,
            reference_pos,
            min_base_quality,
            entries,
            pileups: PyOnceLock::new(),
        }
    }

    /// Whether this is the pileup at a position of a contig.
    pub(crate) fn is_at(&self, py: Python<'_>, contig: &str, position: i64) -> bool {
        self.reference_pos == position
            && self
                .reference_name
                .bind(py)
                .to_cow()
                .is_ok_and(|name| name == contig)
    }

    /// The contig and position of this pileup, as `contig:position`.
    pub(crate) fn locus(&self, py: Python<'_>) -> String {
        format!("{}:{}", self.reference_name.bind(py), self.reference_pos)
    }

    fn pileups_of<'py>(&self, py: Python<'py>) -> PyResult<&Bound<'py, PyTuple>> {
        let pileups = self.pileups.get_or_try_init(py, || {
            let reads = self
                .entries
                .iter()
                .map(|held| {
                    Py::new(
                        py,
                        PileupRead {
                            held: held.clone(),
                            position: Some(self.reference_pos),
                        },
                    )
                })
                .collect::<PyResult<Vec<_>>>()?;
            PyTuple::new(py, reads).map(Bound::unbind)
        })?;
        Ok(pileups.bind(py))
    }

    fn fields<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        (
            self.reference_name.clone_ref(py),
            self.reference_pos,
            self.pileups_of(py)?.clone(),
            self.min_base_quality,
        )
            .into_pyobject(py)
    }

    fn passing(&self) -> impl Iterator<Item = PyResult<(&Held, bool)>> {
        let floor = self.min_base_quality;
        self.entries
            .iter()
            .map(move |held| Ok((held, held.passes(floor)?)))
    }
}

#[pymethods]
impl Pileup {
    #[new]
    #[pyo3(signature = (
        reference_name,
        reference_pos,
        pileups,
        min_base_quality = i64::from(DEFAULT_MIN_BASE_QUALITY),
    ))]
    fn new(
        py: Python<'_>,
        reference_name: Bound<'_, PyString>,
        reference_pos: i64,
        pileups: &Bound<'_, PyAny>,
        min_base_quality: i64,
    ) -> PyResult<Self> {
        let pileups = PyTuple::new(py, pileups.try_iter()?.collect::<PyResult<Vec<_>>>()?)?;
        let entries = pileups
            .iter()
            .map(|read| {
                read.cast::<PileupRead>()
                    .map(|read| read.get().held.clone())
                    .map_err(|_| PyTypeError::new_err("every entry of a Pileup is a PileupRead"))
            })
            .collect::<PyResult<Vec<_>>>()?;
        let pileup = Self::of(
            reference_name.unbind(),
            reference_pos,
            min_base_quality,
            entries,
        );
        let _ = pileup.pileups.set(py, pileups.unbind());
        Ok(pileup)
    }

    /// Build the pileup at one position from reads in any order.
    ///
    /// Unmapped reads and reads on other contigs are ignored.
    ///
    /// Args:
    ///     alignments: the reads to pile up.
    ///     contig: the name of the contig.
    ///     pos: the 0-based position on the contig.
    ///     min_base_quality: the quality floor of the pileup's filtered views.
    #[classmethod]
    #[pyo3(signature = (alignments, contig, pos, min_base_quality = i64::from(DEFAULT_MIN_BASE_QUALITY)))]
    fn from_alignments(
        _cls: &Bound<'_, PyType>,
        alignments: &Bound<'_, PyAny>,
        contig: Bound<'_, PyString>,
        pos: i64,
        min_base_quality: i64,
    ) -> PyResult<Self> {
        let mut entries = Vec::new();
        let mut footprint = Footprint::default();
        for alignment in alignments.try_iter()? {
            let alignment = alignment?;
            let mut record = bam::Record::default();
            bridge::read(&alignment, &mut record)?;
            let flags = record.flags();
            let (Some(Ok(reference_id)), Some(Ok(start))) =
                (record.reference_sequence_id(), record.alignment_start())
            else {
                continue;
            };
            if flags.is_unmapped() {
                continue;
            }
            let start = usize::from(start) as i64 - 1;
            let placed = footprint
                .fill(start, record.cigar().as_bytes(), record.sequence().len())
                .map_err(|source| {
                    to_python(crate::Error::InvalidRecord {
                        name: crate::pileup::record_name(&record),
                        source,
                    })
                })?;
            if !placed || !alignment.getattr("reference_name")?.eq(&contig)? {
                continue;
            }
            let other_end = other_end(
                &record,
                reference_id,
                start,
                footprint.end,
                MateCigar::Unsearched,
            )
            .ok();
            let query_length = footprint.query_length;
            let record = Arc::new(Bridged {
                object: alignment.unbind(),
                record,
            });
            let mut push = |kind, offset, length| {
                let mut held = Held::located(&record, kind, offset, length);
                held.other_end = other_end;
                held.query_length = Some(query_length);
                entries.push(held);
            };
            match footprint.locate(pos) {
                Some(Located::Base(offset)) => push(EntryKind::Base, offset, 0),
                Some(Located::Deletion(next)) => push(EntryKind::Deletion, next.unwrap_or(NONE), 0),
                Some(Located::Skip) => push(EntryKind::Skip, NONE, 0),
                None => {}
            }
            if let Some(insertion) = footprint.insertion_at(pos) {
                push(EntryKind::Insertion, insertion.offset, insertion.length);
            }
        }
        Ok(Self::of(contig.unbind(), pos, min_base_quality, entries))
    }

    /// The name of the contig.
    #[getter]
    fn reference_name(&self, py: Python<'_>) -> Py<PyString> {
        self.reference_name.clone_ref(py)
    }

    /// The 0-based position on the contig.
    #[getter]
    fn reference_pos(&self) -> i64 {
        self.reference_pos
    }

    /// The entries of the reads at this position.
    #[getter]
    fn pileups<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        self.pileups_of(py).cloned()
    }

    /// The base quality below which bases are left out of the filtered views.
    #[getter]
    fn min_base_quality(&self) -> i64 {
        self.min_base_quality
    }

    /// The number of reads with a base, a deletion, or a skip at this position.
    ///
    /// Quality is ignored, and insertion entries are not counted, as in htslib's column depth.
    #[getter]
    fn unfiltered_depth(&self) -> usize {
        self.entries
            .iter()
            .filter(|held| held.kind != EntryKind::Insertion)
            .count()
    }

    /// The number of reads with a base or a deletion at this position at the quality floor.
    ///
    /// A deletion is judged by the quality of the read's next base, as pysam does. A skip has no
    /// quality, so it is never counted, while pysam judges a skip by its next base too.
    #[getter]
    fn filtered_depth(&self) -> PyResult<usize> {
        let mut depth = 0;
        for passing in self.passing() {
            let (held, passes) = passing?;
            depth += usize::from(passes && held.kind != EntryKind::Insertion);
        }
        Ok(depth)
    }

    /// The upper-cased bases at this position at the quality floor, in the order of `pileups`.
    ///
    /// Only base entries at `min_base_quality` are listed, so the list lines up with `qualities`.
    /// pysam's `get_query_sequences()`, by contrast, lists every entry, with an empty string for a
    /// deletion or a skip.
    #[getter]
    fn bases(&self) -> PyResult<Vec<char>> {
        let mut bases = Vec::with_capacity(self.entries.len());
        for passing in self.passing() {
            let (held, passes) = passing?;
            if passes
                && held.kind == EntryKind::Base
                && let Some(base) = held.base()?
            {
                bases.push(char::from(base));
            }
        }
        Ok(bases)
    }

    /// The base qualities of the bases at this position at the quality floor, as in `bases`.
    ///
    /// pysam's `get_query_qualities()`, by contrast, lists every entry, with the quality of the
    /// read's next base for a deletion or a skip.
    #[getter]
    fn qualities<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        let floor = self.min_base_quality;
        let mut qualities = Vec::with_capacity(self.entries.len());
        for held in &self.entries {
            if held.kind == EntryKind::Base
                && let Some(quality) = held.quality()?
                && i64::from(quality) >= floor
            {
                qualities.push(quality);
            }
        }
        PyList::new(py, qualities)
    }

    /// A copy of this pileup with one read per template, by query name.
    ///
    /// The read kept for a template is the first of its name in `pileups` whose entry here is a
    /// base or a deletion at `min_base_quality`, or else the first of its name. A builder fills
    /// `pileups` in input order, so of two passing mates in a coordinate-sorted file, the one
    /// that starts first is kept. A mate's skip, or its base under the floor, therefore never
    /// hides the other mate's base, as in htslib, and as in fgbio, whose floor drops a failing
    /// base before its `withoutOverlaps` keeps the first entry. Every entry of the kept read
    /// stays, its insertion entry included, where fgbio keeps only the first entry of each
    /// template.
    fn without_overlaps(&self, py: Python<'_>) -> PyResult<Self> {
        let passing = self.passing().collect::<PyResult<Vec<_>>>()?;
        let kept = templates_kept(&passing, |(held, passes)| {
            (
                held.template(),
                held.record.object.as_ptr() as usize,
                *passes,
            )
        });
        let entries = self
            .entries
            .iter()
            .zip(&kept)
            .filter(|(_, kept)| **kept)
            .map(|(held, _)| held.clone())
            .collect();
        let pileup = Self::of(
            self.reference_name.clone_ref(py),
            self.reference_pos,
            self.min_base_quality,
            entries,
        );
        if let Some(pileups) = self.pileups.get(py) {
            let reads = pileups
                .bind(py)
                .iter()
                .zip(&kept)
                .filter(|(_, kept)| **kept)
                .map(|(read, _)| read)
                .collect::<Vec<_>>();
            let _ = pileup.pileups.set(py, PyTuple::new(py, reads)?.unbind());
        }
        Ok(pileup)
    }

    fn __eq__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        let Ok(other) = other.cast::<Pileup>() else {
            return Ok(py.NotImplemented());
        };
        let equal = self.fields(py)?.eq(other.get().fields(py)?)?;
        Ok(equal.into_pyobject(py)?.to_owned().into_any().unbind())
    }

    fn __hash__(&self, py: Python<'_>) -> PyResult<isize> {
        self.fields(py)?.hash()
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        Ok(format!(
            "Pileup(reference_name={}, reference_pos={}, pileups={}, min_base_quality={})",
            self.reference_name.bind(py).repr()?,
            self.reference_pos,
            self.pileups_of(py)?.repr()?,
            self.min_base_quality,
        ))
    }
}
