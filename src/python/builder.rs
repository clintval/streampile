//! `StreamingPileupBuilder` over any iterable of pysam `AlignedSegment`s.

use std::collections::HashMap;
use std::io;
use std::num::NonZeroUsize;
use std::sync::{Arc, LazyLock, Mutex, MutexGuard, PoisonError};

use noodles::bam;
use noodles::sam::{
    self,
    alignment::record::Flags,
    header::record::value::{
        Map,
        map::{ReferenceSequence, header::tag::SORT_ORDER},
    },
};
use pyo3::PyTraverseError;
use pyo3::exceptions::{PyAttributeError, PyValueError};
use pyo3::gc::PyVisit;
use pyo3::prelude::*;
use pyo3::types::{PyIterator, PyString};

use super::bridge;
use super::pileup::{Held, Pileup};
use super::{Int, to_python};
use crate::pileup::Derived;
use crate::source::{AlignmentRecord, RecordSource};
use crate::{DEFAULT_EXCLUDE_FLAGS, DEFAULT_MIN_BASE_QUALITY, Error};

static EMPTY: LazyLock<bam::Record> = LazyLock::new(bam::Record::default);

/// A pysam record, the BAM record copied from it when it was read, and what is worked out about
/// the copy when first asked for.
pub(crate) struct Bridged {
    pub object: Py<PyAny>,
    pub record: bam::Record,
    pub derived: Derived,
}

impl Bridged {
    pub(crate) fn new(object: Py<PyAny>, record: bam::Record) -> Self {
        Self {
            object,
            record,
            derived: Derived::default(),
        }
    }
}

/// A record a Python source reads, shared with the pileups that hold it.
#[derive(Default)]
pub(crate) struct PyRecord(Option<Arc<Bridged>>);

impl PyRecord {
    /// The pysam record and its copy, which every record read from a source has.
    pub(crate) fn shared(&self) -> &Arc<Bridged> {
        self.0.as_ref().expect("a record read from a source")
    }

    fn object(&self, py: Python<'_>) -> Py<PyAny> {
        self.shared().object.clone_ref(py)
    }
}

impl AlignmentRecord for PyRecord {
    fn bam(&self) -> &bam::Record {
        self.0.as_ref().map_or(&EMPTY, |shared| &shared.record)
    }

    fn release(&mut self) {
        self.0 = None;
    }
}

/// A Python object held by the builder's class and shared with the source and the callbacks it
/// gives the Rust builder, so that the garbage collector sees the one reference and can clear it.
struct Shared<T = PyAny>(Arc<Mutex<Option<Py<T>>>>);

impl<T> Shared<T> {
    fn new(object: Option<Py<T>>) -> Self {
        Self(Arc::new(Mutex::new(object)))
    }

    fn share(&self) -> Self {
        Self(Arc::clone(&self.0))
    }

    fn get(&self, py: Python<'_>) -> Option<Py<T>> {
        self.lock().as_ref().map(|object| object.clone_ref(py))
    }

    fn take(&self) -> Option<Py<T>> {
        self.lock().take()
    }

    /// Visits the object, unless another thread is reading it, as the garbage collector asks.
    fn visit(&self, visit: &PyVisit<'_>) -> Result<(), PyTraverseError> {
        match self.0.try_lock() {
            Ok(object) => visit.call(object.as_ref()),
            Err(_) => Ok(()),
        }
    }

    fn lock(&self) -> MutexGuard<'_, Option<Py<T>>> {
        self.0.lock().unwrap_or_else(PoisonError::into_inner)
    }
}

/// The records of a Python iterable, and the first one, read already to find the header.
struct PySource {
    records: Shared<PyIterator>,
    first: Shared,
}

impl RecordSource for PySource {
    type Record = PyRecord;

    fn read_record(&mut self, record: &mut PyRecord) -> io::Result<bool> {
        Python::attach(|py| {
            let object = match (self.first.take(), self.records.get(py)) {
                (Some(object), _) => object.into_bound(py),
                (None, Some(records)) => match records.into_bound(py).next() {
                    Some(object) => object.map_err(io::Error::other)?,
                    None => return Ok(false),
                },
                (None, None) => return Ok(false),
            };
            let mut bam = bam::Record::default();
            bridge::read(&object, &mut bam).map_err(io::Error::other)?;
            *record = PyRecord(Some(Arc::new(Bridged::new(object.unbind(), bam))));
            Ok(true)
        })
    }
}

type Builder = crate::StreamingPileupBuilder<'static, PySource>;

/// Build pileups from coordinate-sorted reads in one forward pass.
///
/// Ask for pileups at positions that never move backwards: the same position again returns the
/// pileup already built, and an earlier one raises a `ValueError`. Each read's CIGAR is walked
/// once, when the read is first reached, so building a pileup costs one lookup per read there.
/// Reads are piled up in Rust, from a copy of each read's record made as it is read.
///
/// Every read, filtered or not, is handed to `tap` exactly once and in input order, as soon as
/// the builder has moved past it and every read before it, or when the builder closes. A read
/// can therefore be changed, e.g. tagged, while it is in a pileup and then written by `tap`.
/// Pileups see each read as it was when the builder read it; changes made after that,
/// including in `read_filter`, reach `tap` and `alignment` but not later pileups.
/// Keeping input order means a read is held until every read before it has been passed, so a
/// long read, e.g. one with a long reference skip, holds back every read that starts within it:
/// with a `tap`, buffering behind the longest active read is inherent. Without a `tap`, a read
/// is dropped as soon as the builder has moved past it.
///
/// An exception while advancing to a position, from `records`, a malformed read, `read_filter`,
/// or `tap`, stops the builder, so a pileup never misses a read: every later pileup raises a
/// `ValueError`, and closing still hands every read not yet handed over to `tap`.
///
/// A builder can be used from any thread, one call at a time: a call made while another runs,
/// such as from `read_filter` or `tap`, raises a `RuntimeError`.
///
/// ```python
/// with (
///     AlignmentFile("in.bam", threads=4) as source,
///     AlignmentFile("out.bam", "wb", template=source) as sink,
///     StreamingPileupBuilder(source, tap=sink.write) as builder,
/// ):
///     pileup = builder.pileup("chr1", 100)
/// ```
#[pyclass(module = "streampile", name = "StreamingPileupBuilder")]
pub(crate) struct StreamingPileupBuilder {
    inner: Mutex<Builder>,
    header: Option<Py<PyAny>>,
    names: Vec<Py<PyString>>,
    ids: HashMap<String, usize>,
    min_mapping_quality: u8,
    exclude_flags: u16,
    min_base_quality: u8,
    proper_pairs_only: bool,
    records: Shared<PyIterator>,
    first: Shared,
    read_filter: Shared,
    tap: Shared,
    previous: Option<Py<Pileup>>,
    asked: Option<(usize, i64)>,
    closed: bool,
}

#[pymethods]
impl StreamingPileupBuilder {
    /// Start a builder over coordinate-sorted reads.
    ///
    /// The header is the `header` of `records` when it has one, such as an `AlignmentFile`, or
    /// else the first read's, so an empty `AlignmentFile` is checked too. With no header at all,
    /// e.g. from an empty list, every pileup is empty, and contigs are ordered as they are first
    /// asked for.
    ///
    /// Args:
    ///     records: coordinate-sorted reads, such as an open `AlignmentFile`.
    ///     min_mapping_quality: the lowest mapping quality of a read to pile up.
    ///     exclude_flags: reads with any of these SAM flags are not piled up: by default,
    ///         secondary, QC-fail, duplicate, and supplementary reads, as in `tabulate`. htslib
    ///         keeps supplementary reads and fgbio keeps QC-fail reads.
    ///     min_base_quality: the quality floor of each pileup's filtered views: 13 by default,
    ///         as in pysam's `pileup()` and `samtools mpileup`, where `tabulate` counts bases of
    ///         any quality by default.
    ///     proper_pairs_only: pile up only reads flagged as in a proper pair.
    ///     read_filter: a function that keeps a read for pileups when it returns True, e.g.
    ///         `lambda read: read.is_proper_pair`, asked only of reads that pass the other
    ///         filters. A read it rejects still goes to `tap`.
    ///     tap: a function given every read once the builder has moved past it.
    ///
    /// Raises:
    ///     ValueError: if the header does not declare coordinate order, or a quality is not from 0
    ///         to 255 or `exclude_flags` from 0 to 65535.
    #[new]
    #[pyo3(signature = (
        records,
        *,
        min_mapping_quality = Int::from(0),
        exclude_flags = Int::from(i64::from(DEFAULT_EXCLUDE_FLAGS.bits())),
        min_base_quality = Int::from(i64::from(DEFAULT_MIN_BASE_QUALITY)),
        proper_pairs_only = false,
        read_filter = None,
        tap = None,
    ))]
    #[pyo3(
        text_signature = "(records, *, min_mapping_quality=0, exclude_flags=3840, \
        min_base_quality=13, proper_pairs_only=False, read_filter=None, tap=None)"
    )]
    #[allow(clippy::too_many_arguments)]
    fn new(
        py: Python<'_>,
        records: &Bound<'_, PyAny>,
        min_mapping_quality: Int,
        exclude_flags: Int,
        min_base_quality: Int,
        proper_pairs_only: bool,
        read_filter: Option<Py<PyAny>>,
        tap: Option<Py<PyAny>>,
    ) -> PyResult<Self> {
        let min_mapping_quality = min_mapping_quality.of("min_mapping_quality", u8::MAX)?;
        let exclude_flags = exclude_flags.of("exclude_flags", u16::MAX)?;
        let min_base_quality = min_base_quality.of("min_base_quality", u8::MAX)?;
        let source_header = match records.getattr("header") {
            Ok(header) => Some(header),
            Err(error) if error.is_instance_of::<PyAttributeError>(py) => None,
            Err(error) => return Err(error),
        };
        let iterator = records.try_iter()?;
        let first = iterator.clone().next().transpose()?;
        let alignment_header = py.import("pysam")?.getattr("AlignmentHeader")?;
        let header = match source_header {
            Some(header) if header.is_instance(&alignment_header)? => Some(header),
            _ => match &first {
                Some(first) => Some(first.getattr("header")?).filter(|header| !header.is_none()),
                None => None,
            },
        };
        let mut names = Vec::new();
        if let Some(header) = &header {
            if sort_order(&header.str()?.to_cow()?) != Some("coordinate") {
                return Err(to_python(Error::NotCoordinateSorted { found: None }));
            }
            for name in header.getattr("references")?.try_iter()? {
                names.push(name?.cast_into::<PyString>()?.unbind());
            }
        }
        let lengths: Vec<usize> = match &header {
            Some(header) => header.getattr("lengths")?.extract()?,
            None => Vec::new(),
        };
        let mut sam_header = sam::Header::builder().set_header({
            let mut map = Map::<sam::header::record::value::map::Header>::default();
            map.other_fields_mut()
                .insert(SORT_ORDER, "coordinate".into());
            map
        });
        let mut ids = HashMap::with_capacity(names.len());
        for (id, (name, length)) in names.iter().zip(lengths).enumerate() {
            let name = name.bind(py).to_str()?.to_owned();
            let length = NonZeroUsize::new(length).unwrap_or(NonZeroUsize::MIN);
            sam_header = sam_header
                .add_reference_sequence(name.as_bytes(), Map::<ReferenceSequence>::new(length));
            ids.entry(name).or_insert(id);
        }
        let records = Shared::new(Some(iterator.unbind()));
        let first = Shared::new(first.map(Bound::unbind));
        let source = PySource {
            records: records.share(),
            first: first.share(),
        };
        let read_filter = Shared::new(read_filter);
        let tap = Shared::new(tap);
        let inner = Builder::new(source, &sam_header.build())
            .map_err(to_python)?
            .min_mapping_quality(min_mapping_quality)
            .exclude_flags(Flags::from_bits_retain(exclude_flags))
            .min_base_quality(min_base_quality)
            .proper_pairs_only(proper_pairs_only);
        let inner = with_callbacks(inner, &read_filter, &tap);
        Ok(Self {
            inner: Mutex::new(inner),
            header: header.map(Bound::unbind),
            names,
            ids,
            min_mapping_quality,
            exclude_flags,
            min_base_quality,
            proper_pairs_only,
            records,
            first,
            read_filter,
            tap,
            previous: None,
            asked: None,
            closed: false,
        })
    }

    fn __enter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    #[allow(clippy::needless_pass_by_value)]
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        visit.call(&self.header)?;
        visit.call(&self.previous)?;
        self.records.visit(&visit)?;
        self.first.visit(&visit)?;
        self.read_filter.visit(&visit)?;
        self.tap.visit(&visit)
    }

    fn __clear__(&mut self) {
        exclusive(&mut self.inner).abandon();
        self.header = None;
        self.previous = None;
        self.records.take();
        self.first.take();
        self.read_filter.take();
        self.tap.take();
    }

    #[pyo3(signature = (*_exception))]
    fn __exit__(&mut self, _exception: &Bound<'_, pyo3::types::PyTuple>) -> PyResult<()> {
        self.close()
    }

    /// Stop, first handing every read not yet handed to `tap` to it, in input order.
    ///
    /// With a `tap`, the rest of the input is read to the end, so an output written by `tap` is
    /// complete. Without one, no more of the input is read. A read goes to `tap` once even when
    /// `tap` raises, and closing again hands over the reads after it.
    fn close(&mut self) -> PyResult<()> {
        self.closed = true;
        exclusive(&mut self.inner).close().map_err(to_python)
    }

    /// Whether a read passes the built-in filters, is placed, and then passes `read_filter`.
    fn accepts(&mut self, record: &Bound<'_, PyAny>) -> PyResult<bool> {
        let mut bam = bam::Record::default();
        bridge::read(record, &mut bam)?;
        let record = PyRecord(Some(Arc::new(Bridged::new(record.clone().unbind(), bam))));
        exclusive(&mut self.inner)
            .accepts(&record)
            .map_err(to_python)
    }

    /// Advance to a position, at or after the last one, and pile up the reads there.
    ///
    /// Args:
    ///     contig: the name of the contig.
    ///     pos: the 0-based position on the contig.
    ///
    /// Raises:
    ///     ValueError: if the builder is closed or stopped by an earlier exception, the contig is
    ///         not in the header, the position is negative, or the position is before the last
    ///         one asked for.
    fn pileup(&mut self, py: Python<'_>, contig: &str, pos: i64) -> PyResult<Py<Pileup>> {
        if self.closed {
            return Err(to_python(Error::Closed));
        }
        if let Some(previous) = &self.previous
            && previous.get().is_at(py, contig, pos)
        {
            return Ok(previous.clone_ref(py));
        }
        if pos < 0 {
            return Err(PyValueError::new_err(format!(
                "Position must be non-negative, found: {pos}"
            )));
        }
        let pileup = if self.header.is_some() {
            let id = *self
                .ids
                .get(contig)
                .ok_or_else(|| to_python(Error::UnknownContig(contig.to_owned())))?;
            let pileup = match exclusive(&mut self.inner).pileup_at(id, pos as usize) {
                Ok(pileup) => pileup,
                Err(Error::Backwards { .. }) => return Err(self.backwards(py, contig, pos)),
                Err(error) => {
                    self.previous = None;
                    return Err(to_python(error));
                }
            };
            let entries = pileup.iter().map(|entry| Held::of(&entry, pos)).collect();
            Pileup::of(
                self.names[id].clone_ref(py),
                pos,
                i64::from(self.min_base_quality),
                entries,
            )
        } else {
            let next = self.ids.len();
            let id = *self.ids.entry(contig.to_owned()).or_insert(next);
            if self.asked.is_some_and(|asked| (id, pos) < asked) {
                return Err(self.backwards(py, contig, pos));
            }
            self.asked = Some((id, pos));
            Pileup::of(
                PyString::new(py, contig).unbind(),
                pos,
                i64::from(self.min_base_quality),
                Vec::new(),
            )
        };
        let pileup = Py::new(py, pileup)?;
        self.previous = Some(pileup.clone_ref(py));
        Ok(pileup)
    }

    /// Yield the pileup at every position from `start` to `end`, covered or not.
    ///
    /// Args:
    ///     contig: the name of the contig.
    ///     start: the 0-based first position.
    ///     end: the 0-based position after the last one.
    ///
    /// Raises:
    ///     ValueError: if `end` is before `start`, or as `pileup` does.
    fn columns(slf: Bound<'_, Self>, contig: String, start: i64, end: i64) -> Columns {
        Columns {
            builder: slf.unbind(),
            contig,
            start,
            next: start,
            end,
        }
    }

    /// The header of the reads, or `None` without one.
    #[getter]
    fn header(&self, py: Python<'_>) -> Option<Py<PyAny>> {
        self.header.as_ref().map(|header| header.clone_ref(py))
    }

    /// The lowest mapping quality of a read to pile up.
    #[getter]
    fn min_mapping_quality(&self) -> u8 {
        self.min_mapping_quality
    }

    /// The SAM flags of reads that are not piled up.
    #[getter]
    fn exclude_flags(&self) -> u16 {
        self.exclude_flags
    }

    /// The quality floor of each pileup's filtered views.
    #[getter]
    fn min_base_quality(&self) -> u8 {
        self.min_base_quality
    }

    /// Whether only reads flagged as in a proper pair are piled up.
    #[getter]
    fn proper_pairs_only(&self) -> bool {
        self.proper_pairs_only
    }

    /// The function that keeps a read for pileups, if any.
    #[getter]
    fn read_filter(&self, py: Python<'_>) -> Option<Py<PyAny>> {
        self.read_filter.get(py)
    }

    /// The function given every read once the builder has moved past it, if any.
    #[getter]
    fn tap(&self, py: Python<'_>) -> Option<Py<PyAny>> {
        self.tap.get(py)
    }

    /// The last pileup built, which a repeated position returns again.
    #[getter]
    fn previous_pileup(&self, py: Python<'_>) -> Option<Py<Pileup>> {
        self.previous.as_ref().map(|pileup| pileup.clone_ref(py))
    }
}

impl Drop for StreamingPileupBuilder {
    fn drop(&mut self) {
        exclusive(&mut self.inner).abandon();
    }
}

impl StreamingPileupBuilder {
    fn backwards(&self, py: Python<'_>, contig: &str, pos: i64) -> PyErr {
        let from = self
            .previous
            .as_ref()
            .map_or_else(String::new, |previous| previous.get().locus(py));
        PyValueError::new_err(format!(
            "Attempted to advance to {contig}:{pos} from {from}."
        ))
    }
}

/// The pileup at every position of a span, from `StreamingPileupBuilder.columns`.
#[pyclass(module = "streampile")]
pub(crate) struct Columns {
    builder: Py<StreamingPileupBuilder>,
    contig: String,
    start: i64,
    next: i64,
    end: i64,
}

#[pymethods]
impl Columns {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<Py<Pileup>>> {
        if self.end < self.start {
            return Err(PyValueError::new_err(format!(
                "End {} is before start {}.",
                self.end, self.start
            )));
        }
        if self.next >= self.end {
            return Ok(None);
        }
        let position = self.next;
        self.next += 1;
        self.builder
            .bind(py)
            .try_borrow_mut()?
            .pileup(py, &self.contig, position)
            .map(Some)
    }
}

/// A Rust builder that asks the Python `read_filter` and hands records to the Python `tap`, where
/// they are given.
fn with_callbacks(mut inner: Builder, read_filter: &Shared, tap: &Shared) -> Builder {
    if read_filter.lock().is_some() {
        let read_filter = read_filter.share();
        inner = inner.try_read_filter(move |record: &PyRecord| {
            Python::attach(|py| match read_filter.get(py) {
                Some(read_filter) => read_filter
                    .bind(py)
                    .call1((record.object(py),))
                    .and_then(|kept| kept.is_truthy())
                    .map_err(io::Error::other),
                None => Ok(true),
            })
        });
    }
    if tap.lock().is_some() {
        let tap = tap.share();
        inner = inner.tap(move |record: PyRecord| {
            Python::attach(|py| match tap.get(py) {
                Some(tap) => tap
                    .bind(py)
                    .call1((record.object(py),))
                    .map(drop)
                    .map_err(io::Error::other),
                None => Ok(()),
            })
        });
    }
    inner
}

/// The Rust builder, borrowed without locking: PyO3 lends a builder mutably to one caller at a
/// time, and the mutex only makes the class shareable between threads.
fn exclusive(inner: &mut Mutex<Builder>) -> &mut Builder {
    inner.get_mut().unwrap_or_else(PoisonError::into_inner)
}

/// The sort order an `@HD` line of a SAM header declares, if any.
fn sort_order(text: &str) -> Option<&str> {
    let line = text.lines().find(|line| line.starts_with("@HD\t"))?;
    line.split('\t')
        .filter_map(|field| field.strip_prefix("SO:"))
        .next_back()
}
