//! `Pileup`, `PileupRead`, and `PileupTemplate`, owned snapshots of a column that outlive the
//! builder's next move.

use std::sync::Arc;

use noodles::bam;
use noodles::sam::alignment::record::Cigar as _;
use noodles::sam::alignment::record::cigar::op::Kind;
use pyo3::exceptions::{PyIndexError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::pyclass::CompareOp;
use pyo3::sync::PyOnceLock;
use pyo3::types::{IntoPyDict, PyDict, PyList, PyString, PyTuple, PyType};

use super::builder::Bridged;
use super::{Int, bridge, to_python};
use noodles::sam::alignment::record::Flags;

use crate::footprint::Footprint;
use crate::overlap::{AgreementStrategy, DisagreementStrategy, Observation};
use crate::pileup::{
    EntryKind, MISSING_BASE_QUALITY, NONE, Template, TemplateName, TemplateRead, name_of,
    number_templates, templates,
};
use crate::template::{entry_from_five_prime, from_five_prime};
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

const PILEUP_FIELDS: [&str; 4] = [
    "reference_name",
    "reference_pos",
    "pileups",
    "min_base_quality",
];

static PILEUP_READ_TYPES: PyOnceLock<[Py<PyAny>; 4]> = PyOnceLock::new();

/// One read at one position: the record, the fields of a `PileupRead`, the position where a
/// pileup made it, and what is worked out about the record once, where a builder made it.
#[derive(Clone)]
pub(crate) struct Held {
    record: Arc<Bridged>,
    kind: EntryKind,
    query_position: i64,
    query_position_or_next: i64,
    insertion_offset: i64,
    insertion_length: i64,
    position: Option<i64>,
    five_prime_distance: Option<usize>,
}

impl Held {
    /// The entry of a builder's pileup at a position.
    pub(crate) fn of(entry: &PileupEntry<'_, super::builder::PyRecord>, position: i64) -> Self {
        let present = |value: Option<usize>| value.map_or(ABSENT, |value| value as i64);
        Self {
            record: Arc::clone(entry.source_record().shared()),
            kind: entry.kind(),
            query_position: present(entry.query_position()),
            query_position_or_next: present(entry.query_position_or_next()),
            insertion_offset: present(entry.insertion_offset()),
            insertion_length: entry.insertion_length() as i64,
            position: Some(position),
            five_prime_distance: entry.five_prime_distance(),
        }
    }

    /// An entry of a read of a query length placed by its footprint at a position, as a builder
    /// would make it.
    fn located(
        record: &Arc<Bridged>,
        query_length: u32,
        position: i64,
        (kind, offset, length): (EntryKind, u32, u32),
    ) -> Self {
        let five_prime_distance = entry_from_five_prime(
            kind,
            record.record.flags().is_reverse_complemented(),
            query_length as usize,
            (offset != NONE).then_some(offset as usize),
        );
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
            position: Some(position),
            five_prime_distance,
        }
    }

    /// An entry made by hand from the fields of a `PileupRead`.
    fn by_hand(
        record: Arc<Bridged>,
        kind: EntryKind,
        [
            query_position,
            query_position_or_next,
            insertion_offset,
            insertion_length,
        ]: [i64; 4],
    ) -> Self {
        let bam = &record.record;
        let five_prime_distance = bam
            .cigar()
            .read_length()
            .ok()
            .filter(|_| kind == EntryKind::Base)
            .and_then(|length| {
                let reverse = bam.flags().is_reverse_complemented();
                from_five_prime(reverse, length, usize::try_from(query_position).ok()?)
            });
        Self {
            record,
            kind,
            query_position,
            query_position_or_next,
            insertion_offset,
            insertion_length,
            position: None,
            five_prime_distance,
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

    fn template_end_distance(&self) -> crate::Result<Option<usize>> {
        let position = self.position.or_else(|| self.position_of_base());
        let Some(position) = position.and_then(|position| usize::try_from(position).ok()) else {
            return Ok(None);
        };
        let ends = self.record.derived.ends(self.bam())?;
        Ok(ends.and_then(|ends| ends.distance(position)))
    }

    fn is_fr_pair(&self) -> crate::Result<bool> {
        self.record.derived.is_fr_pair(self.bam())
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

    fn template(&self) -> TemplateName<'_> {
        TemplateName::of(
            self.bam(),
            &self.record.derived,
            Arc::as_ptr(&self.record) as usize,
        )
    }

    /// What the read holds here, as its template's vote counts it.
    fn observation(&self) -> PyResult<Observation> {
        Ok(Observation {
            kind: self.kind,
            base: self.base()?,
            quality: self.quality()?,
        })
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

impl TemplateRead for Held {
    fn flags(&self) -> Flags {
        self.bam().flags()
    }

    fn five_prime_distance(&self) -> Option<usize> {
        self.five_prime_distance
    }

    fn template_end_distance(&self) -> crate::Result<Option<usize>> {
        Held::template_end_distance(self)
    }

    fn is_fr_pair(&self) -> crate::Result<bool> {
        Held::is_fr_pair(self)
    }
}

/// The fields of a `PileupRead`, or a tuple, or `None` for anything else.
fn as_tuple<'py>(
    py: Python<'py>,
    value: &Bound<'py, PyAny>,
) -> PyResult<Option<Bound<'py, PyTuple>>> {
    if let Ok(read) = value.cast::<PileupRead>() {
        return read.get().held.fields(py).map(Some);
    }
    Ok(value.cast::<PyTuple>().ok().cloned())
}

/// A pysam record and the BAM record copied from it.
fn bridged(alignment: &Bound<'_, PyAny>) -> PyResult<Arc<Bridged>> {
    let mut record = bam::Record::default();
    bridge::read(alignment, &mut record)?;
    Ok(Arc::new(Bridged::new(alignment.clone().unbind(), record)))
}

/// A tuple of the `PileupRead`s of entries.
fn reads_of<'a>(
    py: Python<'_>,
    entries: impl IntoIterator<Item = &'a Held>,
) -> PyResult<Py<PyTuple>> {
    let reads = entries
        .into_iter()
        .map(|held| Py::new(py, PileupRead { held: held.clone() }))
        .collect::<PyResult<Vec<_>>>()?;
    Ok(PyTuple::new(py, reads)?.unbind())
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
/// A `PileupRead` behaves as the `NamedTuple` of its six fields, though it is not a `tuple`: it
/// unpacks, indexes, compares, hashes, matches, copies, and pickles as one, and has `_make`,
/// `_asdict`, and `_replace`. Its bases and qualities are those of the read as it was when it was
/// piled up, while `alignment` is the very object given to the builder, so it can be changed,
/// e.g. tagged, and written on.
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
        Ok(Self {
            held: Held::by_hand(
                bridged(alignment)?,
                kind_of(pileup_type)?,
                [
                    optional(query_position),
                    optional(query_position_or_next),
                    optional(insertion_offset),
                    insertion_length,
                ],
            ),
        })
    }

    #[classattr]
    #[pyo3(name = "_fields")]
    fn field_names(py: Python<'_>) -> PyResult<Bound<'_, PyTuple>> {
        PyTuple::new(py, FIELDS)
    }

    #[classattr]
    fn __match_args__(py: Python<'_>) -> PyResult<Bound<'_, PyTuple>> {
        PyTuple::new(py, FIELDS)
    }

    #[classattr]
    #[pyo3(name = "_field_defaults")]
    fn field_defaults(py: Python<'_>) -> PyResult<Bound<'_, PyDict>> {
        let defaults = PyDict::new(py);
        defaults.set_item("insertion_offset", py.None())?;
        defaults.set_item("insertion_length", 0)?;
        Ok(defaults)
    }

    #[classmethod]
    #[pyo3(name = "_make")]
    fn make<'py>(
        cls: &Bound<'py, PyType>,
        iterable: &Bound<'py, PyAny>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let fields = iterable.try_iter()?.collect::<PyResult<Vec<_>>>()?;
        cls.call1(PyTuple::new(cls.py(), fields)?)
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
    /// `positionInReadInReadOrder` minus one. For a deletion or a skip, which holds no base, it is
    /// the number of the read's bases sequenced before the position, as `template_end_distance`
    /// counts the bases after it. It is `None` for an insertion entry, and for an entry made by hand
    /// that holds no base, whose position is unknown.
    #[getter]
    fn five_prime_distance(&self) -> Option<usize> {
        self.held.five_prime_distance()
    }

    /// The number of the template's bases between the position and the template's other end, the
    /// 5′ end of the mate of a read in an FR pair: 0 at the mate's 5′ end.
    ///
    /// It walks the read's CIGAR and the mate's, from its `MC` tag, so an indel counts by its
    /// length; soft-clipped bases count and hard-clipped bases, absent from the records, do not.
    /// Where the mate has no base, a position both reads align to carries the count, or else the
    /// reference between them; the template length (TLEN) is never read. It is `None` for a read
    /// that `is_fr_pair` does not call a read of an FR pair, for a read of an FR pair only at a
    /// position past the mate's 5′ end, where a read runs through its mate, and for an entry made
    /// by hand that holds no base, whose position is unknown.
    ///
    /// Raises:
    ///     ValueError: for a read of an FR pair with no `MC` tag, or one that is not a CIGAR
    ///         string spanning at least one base.
    #[getter]
    fn template_end_distance(&self) -> PyResult<Option<usize>> {
        self.held.template_end_distance().map_err(to_python)
    }

    /// Whether the read is a read of an FR pair, as htsjdk 5.0.0's `getPairOrientation` says.
    ///
    /// A pair is FR when its reads are paired, mapped to one contig, and on opposite strands, and
    /// the forward read's aligned 5′ position is at or before the reverse read's, so a pair whose
    /// 5′ ends coincide is FR. A forward read takes its mate's aligned end from its `MC` tag and
    /// never from the template length (TLEN).
    ///
    /// Raises:
    ///     ValueError: for a forward read of a pair otherwise FR with no `MC` tag, or one that is
    ///         not a CIGAR string.
    #[getter]
    fn is_fr_pair(&self) -> PyResult<bool> {
        self.held.is_fr_pair().map_err(to_python)
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
        let mut moved = false;
        if let Some(changes) = changes {
            for (name, value) in changes {
                if !FIELDS.contains(&name.extract::<&str>()?) {
                    return Err(PyValueError::new_err(format!(
                        "Got unexpected field names: [{}]",
                        name.repr()?
                    )));
                }
                moved |= name.eq("alignment")?;
                fields.set_item(name, value)?;
            }
        }
        let field = |name: &str| fields.as_any().get_item(name);
        let record = if moved {
            bridged(&field("alignment")?)?
        } else {
            Arc::clone(&self.held.record)
        };
        let mut held = Held::by_hand(
            record,
            kind_of(&field("pileup_type")?.extract::<String>()?)?,
            [
                optional(field("query_position")?.extract()?),
                optional(field("query_position_or_next")?.extract()?),
                optional(field("insertion_offset")?.extract()?),
                field("insertion_length")?.extract()?,
            ],
        );
        if !moved {
            held.position = self.held.position;
            let placed =
                |held: &Held| (held.kind, held.query_position, held.query_position_or_next);
            if placed(&held) == placed(&self.held) {
                held.five_prime_distance = self.held.five_prime_distance;
            }
        }
        Ok(Self { held })
    }

    #[pyo3(signature = (**changes))]
    fn __replace__(&self, py: Python<'_>, changes: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
        self.replace(py, changes)
    }

    fn count<'py>(
        &self,
        py: Python<'py>,
        value: &Bound<'py, PyAny>,
    ) -> PyResult<Bound<'py, PyAny>> {
        self.held.fields(py)?.call_method1("count", (value,))
    }

    #[pyo3(signature = (value, *bounds))]
    fn index<'py>(
        &self,
        py: Python<'py>,
        value: &Bound<'py, PyAny>,
        bounds: &Bound<'py, PyTuple>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let arguments: Vec<_> = std::iter::once(value.clone()).chain(bounds).collect();
        let arguments = PyTuple::new(py, arguments)?;
        self.held.fields(py)?.call_method1("index", arguments)
    }

    fn __contains__(&self, py: Python<'_>, value: &Bound<'_, PyAny>) -> PyResult<bool> {
        self.held.fields(py)?.contains(value)
    }

    fn __add__<'py>(&self, py: Python<'py>, other: &Bound<'py, PyAny>) -> PyResult<Py<PyAny>> {
        let Some(other) = as_tuple(py, other)? else {
            return Ok(py.NotImplemented());
        };
        Ok(self.held.fields(py)?.add(other)?.unbind())
    }

    fn __radd__<'py>(&self, py: Python<'py>, other: &Bound<'py, PyAny>) -> PyResult<Py<PyAny>> {
        let Some(other) = as_tuple(py, other)? else {
            return Ok(py.NotImplemented());
        };
        Ok(other.add(self.held.fields(py)?)?.unbind())
    }

    fn __mul__<'py>(&self, py: Python<'py>, times: &Bound<'py, PyAny>) -> PyResult<Py<PyAny>> {
        match self.held.fields(py)?.as_any().mul(times) {
            Ok(repeated) => Ok(repeated.unbind()),
            Err(_) => Ok(py.NotImplemented()),
        }
    }

    fn __rmul__<'py>(&self, py: Python<'py>, times: &Bound<'py, PyAny>) -> PyResult<Py<PyAny>> {
        self.__mul__(py, times)
    }

    fn __getnewargs__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        self.held.fields(py)
    }

    fn __reduce__<'py>(
        slf: &Bound<'py, Self>,
    ) -> PyResult<(Bound<'py, PyType>, Bound<'py, PyTuple>)> {
        Ok((slf.get_type(), slf.get().held.fields(slf.py())?))
    }

    fn __copy__(slf: Bound<'_, Self>) -> Bound<'_, Self> {
        slf
    }

    fn __deepcopy__(&self, py: Python<'_>, memo: &Bound<'_, PyAny>) -> PyResult<Self> {
        let alignment = py
            .import("copy")?
            .call_method1("deepcopy", (self.held.record.object.bind(py), memo))?;
        let record = Bridged::new(alignment.unbind(), self.held.record.record.clone());
        Ok(Self {
            held: Held {
                record: Arc::new(record),
                ..self.held.clone()
            },
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
        let Some(other) = as_tuple(py, other)? else {
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
/// worked out from it in Rust. It behaves as a frozen dataclass of its four fields: it compares,
/// hashes, matches, copies, and pickles by them, and `dataclasses.replace`, `fields`, `asdict`,
/// and `astuple` take it.
///
/// Attributes:
///     reference_name: the name of the contig.
///     reference_pos: the 0-based position on the contig.
///     pileups: the entries of the reads at this position.
///     min_base_quality: the base quality below which bases are left out of `filtered_depth`,
///         `bases`, and `qualities`.
#[pyclass(module = "streampile", name = "Pileup", frozen, weakref)]
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
        let pileups = self
            .pileups
            .get_or_try_init(py, || reads_of(py, &self.entries))?;
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
        min_base_quality = Int::from(i64::from(DEFAULT_MIN_BASE_QUALITY)),
    ))]
    #[pyo3(text_signature = "(reference_name, reference_pos, pileups, min_base_quality=13)")]
    fn new(
        py: Python<'_>,
        reference_name: Bound<'_, PyString>,
        reference_pos: i64,
        pileups: &Bound<'_, PyAny>,
        min_base_quality: Int,
    ) -> PyResult<Self> {
        let min_base_quality = i64::from(min_base_quality.of("min_base_quality", u8::MAX)?);
        let pileups = PyTuple::new(py, pileups.try_iter()?.collect::<PyResult<Vec<_>>>()?)?;
        let entries = pileups
            .iter()
            .map(|read| {
                read.cast::<PileupRead>()
                    .map(|read| Held {
                        position: Some(reference_pos),
                        ..read.get().held.clone()
                    })
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
    #[pyo3(signature = (
        alignments,
        contig,
        pos,
        min_base_quality = Int::from(i64::from(DEFAULT_MIN_BASE_QUALITY)),
    ))]
    #[pyo3(text_signature = "(alignments, contig, pos, min_base_quality=13)")]
    fn from_alignments(
        _cls: &Bound<'_, PyType>,
        alignments: &Bound<'_, PyAny>,
        contig: Bound<'_, PyString>,
        pos: i64,
        min_base_quality: Int,
    ) -> PyResult<Self> {
        let min_base_quality = i64::from(min_base_quality.of("min_base_quality", u8::MAX)?);
        let mut entries = Vec::new();
        let mut footprint = Footprint::default();
        for alignment in alignments.try_iter()? {
            let alignment = alignment?;
            if !alignment.getattr("reference_name")?.eq(&contig)? {
                continue;
            }
            let mut record = bam::Record::default();
            bridge::read(&alignment, &mut record)?;
            let flags = record.flags();
            let (Some(Ok(_)), Some(Ok(start))) =
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
            if !placed {
                continue;
            }
            let query_length = footprint.query_length;
            let record = Arc::new(Bridged::new(alignment.unbind(), record));
            footprint.entries_at(pos, |kind, offset, length| {
                entries.push(Held::located(
                    &record,
                    query_length,
                    pos,
                    (kind, offset, length),
                ));
            });
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

    /// One observation per template at this position, its reads grouped by query name, in the
    /// order of each template's first entry in `pileups`. A read with no name is a template of its
    /// own.
    ///
    /// Where two reads of a template hold bases, they are called into one as fgbio's
    /// `CallOverlappingConsensusBases` calls them, as fgumi implements it: `agreement` makes the
    /// quality of equal bases and `disagreement` the base and quality of different ones. Where the
    /// strategies leave a read's base or quality unchanged, the template takes the higher quality.
    ///
    /// fgumi defines no more than that, so at each position a no-call (`N`) is left alone, as fgumi
    /// leaves it: the other read's base stands at its own quality, and two no-calls are an `N` at
    /// the higher quality. A no-call is still a base, so it stands over the other read's deletion,
    /// which fgumi leaves alone too. A read with a deletion or a skip holds no base, so a template
    /// whose other read holds one has that base at its own quality; a template with no base is a
    /// deletion if either read holds one, at the higher of their qualities, or else a skip.
    /// Insertion entries are no part of a template, so a read whose only entry here is an insertion
    /// adds nothing to its template. A read under `min_base_quality` does not vote, as `bases`,
    /// `qualities`, and `filtered_depth` leave it out: a base, or a deletion judged by its next
    /// base, under the floor. So a mate under the floor neither masks nor lowers the other mate's
    /// base, and a template none of whose reads votes has no `base` or `qual`. A template with more
    /// than two reads here, as when supplementary records are piled up, calls them in the order of
    /// `pileups`.
    ///
    /// Args:
    ///     agreement: how the quality of two reads holding the same base is made.
    ///     disagreement: how the base and quality of two reads holding different bases are made.
    #[pyo3(signature = (*, agreement = "consensus", disagreement = "consensus"))]
    fn templates(&self, agreement: &str, disagreement: &str) -> PyResult<Vec<PileupTemplate>> {
        let agreement: AgreementStrategy = agreement.parse().map_err(to_python)?;
        let disagreement: DisagreementStrategy = disagreement.parse().map_err(to_python)?;
        let numbers = number_templates(&self.entries, |held| {
            (held.kind != EntryKind::Insertion).then(|| held.template())
        });
        let mut entries = Vec::with_capacity(self.entries.len());
        for (held, number) in self.entries.iter().zip(numbers) {
            if let Some(number) = number {
                entries.push((number, held.clone(), held.observation()?));
            }
        }
        let floor = self.min_base_quality;
        Ok(templates(entries, agreement, disagreement, floor)
            .into_iter()
            .map(|template| PileupTemplate {
                template,
                pileups: PyOnceLock::new(),
            })
            .collect())
    }

    #[classattr]
    fn __match_args__(py: Python<'_>) -> PyResult<Bound<'_, PyTuple>> {
        PyTuple::new(py, PILEUP_FIELDS)
    }

    /// The fields of the dataclass a `Pileup` once was, so that `dataclasses.replace`, `fields`,
    /// `asdict`, and `astuple` take a `Pileup`.
    #[classattr]
    fn __dataclass_fields__(py: Python<'_>) -> PyResult<Bound<'_, PyAny>> {
        let dataclasses = py.import("dataclasses")?;
        let builtins = py.import("builtins")?;
        let floor = dataclasses.getattr("field")?.call(
            (),
            Some(&[("default", DEFAULT_MIN_BASE_QUALITY)].into_py_dict(py)?),
        )?;
        let fields = PyList::new(
            py,
            [
                (PILEUP_FIELDS[0], builtins.getattr("str")?)
                    .into_pyobject(py)?
                    .into_any(),
                (PILEUP_FIELDS[1], builtins.getattr("int")?)
                    .into_pyobject(py)?
                    .into_any(),
                (PILEUP_FIELDS[2], builtins.getattr("tuple")?)
                    .into_pyobject(py)?
                    .into_any(),
                (PILEUP_FIELDS[3], builtins.getattr("int")?, floor)
                    .into_pyobject(py)?
                    .into_any(),
            ],
        )?;
        dataclasses
            .call_method1("make_dataclass", ("Pileup", fields))?
            .getattr("__dataclass_fields__")
    }

    #[pyo3(signature = (**changes))]
    fn __replace__<'py>(
        slf: &Bound<'py, Self>,
        changes: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let py = slf.py();
        let fields = PyDict::new(py);
        for (name, value) in PILEUP_FIELDS.iter().zip(slf.get().fields(py)?.iter()) {
            fields.set_item(name, value)?;
        }
        if let Some(changes) = changes {
            fields.update(changes.as_mapping())?;
        }
        slf.get_type().call((), Some(&fields))
    }

    fn __reduce__<'py>(
        slf: &Bound<'py, Self>,
    ) -> PyResult<(Bound<'py, PyType>, Bound<'py, PyTuple>)> {
        Ok((slf.get_type(), slf.get().fields(slf.py())?))
    }

    fn __copy__(slf: Bound<'_, Self>) -> Bound<'_, Self> {
        slf
    }

    fn __deepcopy__(&self, py: Python<'_>, memo: &Bound<'_, PyAny>) -> PyResult<Self> {
        let pileups = py
            .import("copy")?
            .call_method1("deepcopy", (self.pileups_of(py)?, memo))?;
        Self::new(
            py,
            self.reference_name.bind(py).clone(),
            self.reference_pos,
            &pileups,
            Int::from(self.min_base_quality),
        )
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

/// One template at one pileup position: the reads of one query name, their bases called into one.
///
/// A template's strand and distances are those of its first read, the first of a pair or a
/// fragment's only read, worked out from its second read where the first holds no base here.
///
/// Attributes:
///     query_name: the name of the template's reads, `*` for a read with none, which is a
///         template of its own.
///     reads: the entries of the template's reads at the position, usually one or two.
///     pileup_type: whether the template holds a base, a deletion, or a skip.
///     base: the template's base, its reads' bases called into one.
///     qual: the quality of the template's base, or of the next base for a deletion.
#[pyclass(module = "streampile", name = "PileupTemplate", frozen)]
pub(crate) struct PileupTemplate {
    template: Template<Held>,
    pileups: PyOnceLock<Py<PyTuple>>,
}

#[pymethods]
impl PileupTemplate {
    /// The name of the template's reads, `*` for a read with none, which is a template of its own.
    #[getter]
    fn query_name(&self) -> String {
        String::from_utf8_lossy(name_of(self.template.first().bam())).into_owned()
    }

    /// The entries of the template's reads at the position, usually one or two, as in `pileups`.
    #[getter]
    fn reads<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let reads = self
            .pileups
            .get_or_try_init(py, || reads_of(py, self.template.entries()))?;
        Ok(reads.bind(py).clone())
    }

    /// Whether the template holds a base, a deletion, or a skip: a base if a read's base votes, or
    /// else a deletion if a read's deletion votes, or else, with no vote, what its reads hold.
    #[getter]
    fn pileup_type<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        pileup_read_type(py, self.template.called().kind)
    }

    /// The template's upper-cased base, its voting reads' bases called into one, or `None`.
    #[getter]
    fn base(&self) -> Option<char> {
        self.template.called().base.map(char::from)
    }

    /// The quality of the template's base, or of the next base for a deletion, or `None` without
    /// a vote.
    #[getter]
    fn qual(&self) -> Option<u8> {
        self.template.called().quality
    }

    /// Whether the template holds a deletion at the position, and no base that votes.
    #[getter]
    fn is_del(&self) -> bool {
        self.template.called().kind == EntryKind::Deletion
    }

    /// Whether every read of the template skips over the position.
    #[getter]
    fn is_refskip(&self) -> bool {
        self.template.called().kind == EntryKind::Skip
    }

    /// Whether the template's base is a no-call, `N`.
    #[getter]
    fn is_no_call(&self) -> bool {
        self.template.called().base == Some(b'N')
    }

    /// Whether the template's first read is aligned to the reverse strand.
    ///
    /// It is `False` for an F1R2 pair and `True` for an F2R1 pair, read from the second read's
    /// flags without the first.
    #[getter]
    fn is_reverse(&self) -> bool {
        self.template.is_reverse()
    }

    /// The number of the template's bases between its first read's 5′ end and the position.
    ///
    /// It is the first read's `five_prime_distance` where it holds a base here, and otherwise the
    /// second read's `template_end_distance`.
    ///
    /// Raises:
    ///     ValueError: for a read of an FR pair with no usable `MC` tag.
    #[getter]
    fn five_prime_distance(&self) -> PyResult<Option<usize>> {
        self.template.five_prime_distance().map_err(to_python)
    }

    /// The number of the template's bases between the position and its other end, the 5′ end of
    /// the second read of an FR pair.
    ///
    /// It is the first read's `template_end_distance`, and without the first read here the second
    /// read's `five_prime_distance` where that read `is_fr_pair`; it is `None` for any pair that is
    /// not FR.
    ///
    /// Raises:
    ///     ValueError: for a read of an FR pair with no usable `MC` tag.
    #[getter]
    fn template_end_distance(&self) -> PyResult<Option<usize>> {
        self.template.template_end_distance().map_err(to_python)
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        Ok(format!(
            "PileupTemplate(query_name={}, pileup_type={}, base={}, qual={}, reads={})",
            self.query_name().into_pyobject(py)?.repr()?,
            self.pileup_type(py)?.repr()?,
            self.base().into_pyobject(py)?.repr()?,
            self.qual().into_pyobject(py)?.repr()?,
            self.reads(py)?.repr()?,
        ))
    }
}
