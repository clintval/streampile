//! The work of `bin_templates`: the templates that start in each window of a contig.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::io;

use noodles::bam;
use noodles::sam::alignment::record::Flags;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyString;

use super::{Int, bridge, to_python};
use crate::builder::passes;
use crate::pileup::record_name;
use crate::{DEFAULT_EXCLUDE_FLAGS, Error};

/// Secondary and supplementary records, which are never the read a template is counted at.
const NOT_PRIMARY: Flags = Flags::SECONDARY.union(Flags::SUPPLEMENTARY);

/// A contig's place in the order contigs are binned in, and a 0-based position on it.
type Locus = (usize, usize);

/// The templates that start in each window of contigs, binned one contig at a time, in order.
#[pyclass(module = "streampile._native")]
pub(crate) struct TemplateBins {
    alignments: Py<PyAny>,
    contigs: Vec<(Py<PyString>, usize)>,
    ranks: HashMap<usize, usize>,
    window: usize,
    min_mapping_quality: u8,
    exclude_flags: Flags,
    next: usize,
    awaited: BTreeMap<Locus, HashSet<Vec<u8>>>,
    record: bam::Record,
}

#[pymethods]
impl TemplateBins {
    #[new]
    #[pyo3(signature = (
        alignments,
        contigs,
        window,
        *,
        min_mapping_quality = Int::from(0),
        exclude_flags = Int::from(i64::from(DEFAULT_EXCLUDE_FLAGS.bits())),
    ))]
    #[pyo3(
        text_signature = "(alignments, contigs, window, *, min_mapping_quality=0, \
        exclude_flags=3840)"
    )]
    fn new(
        alignments: &Bound<'_, PyAny>,
        contigs: Vec<(Py<PyString>, usize, usize)>,
        window: Int,
        min_mapping_quality: Int,
        exclude_flags: Int,
    ) -> PyResult<Self> {
        let window = window.at_least(1, "window", u32::MAX)? as usize;
        let min_mapping_quality = min_mapping_quality.of("min_mapping_quality", u8::MAX)?;
        let exclude_flags = Flags::from_bits_retain(exclude_flags.of("exclude_flags", u16::MAX)?);
        if !exclude_flags.contains(NOT_PRIMARY) {
            return Err(PyValueError::new_err(format!(
                "exclude_flags must hold the secondary and supplementary flags, {:#x}, found: {:#x}",
                NOT_PRIMARY.bits(),
                exclude_flags.bits()
            )));
        }
        Ok(Self {
            alignments: alignments.clone().unbind(),
            ranks: contigs
                .iter()
                .enumerate()
                .map(|(rank, &(_, id, _))| (id, rank))
                .collect(),
            contigs: contigs
                .into_iter()
                .map(|(name, _, length)| (name, length))
                .collect(),
            window,
            min_mapping_quality,
            exclude_flags,
            next: 0,
            awaited: BTreeMap::new(),
            record: bam::Record::default(),
        })
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<(Py<PyString>, Vec<u64>)>> {
        let Some((contig, length)) = self.contigs.get(self.next) else {
            return Ok(None);
        };
        let (contig, length) = (contig.clone_ref(py), *length);
        let rank = self.next;
        self.next += 1;
        let mut templates = vec![0; length.div_ceil(self.window)];
        let records = self
            .alignments
            .bind(py)
            .call_method1("fetch", (contig.bind(py),))?;
        let mut last = 0;
        for object in records.try_iter()? {
            bridge::read(&object?, &mut self.record)?;
            let record = &self.record;
            let Some(start) = start_of(record.alignment_start())? else {
                continue;
            };
            if start < last {
                return Err(to_python(Error::OutOfOrder {
                    name: record_name(record),
                }));
            }
            last = start;
            let mapping_quality = record.mapping_quality().map_or(255, u8::from);
            if !passes(
                record.flags(),
                mapping_quality,
                self.exclude_flags,
                self.min_mapping_quality,
            ) {
                continue;
            }
            if start >= length {
                return Err(PyValueError::new_err(format!(
                    "Read {} starts past the end of {contig}, which has {length} bases.",
                    record_name(record)
                )));
            }
            let here = (rank, start);
            while let Some(entry) = self.awaited.first_entry()
                && *entry.key() < here
            {
                entry.remove();
            }
            let name = record.name().map(|name| name.to_vec());
            if let Some(name) = &name
                && self
                    .awaited
                    .get_mut(&here)
                    .is_some_and(|names| names.remove(name))
            {
                continue;
            }
            templates[start / self.window] += 1;
            if let Some(name) = name
                && let Some(mate) = mate_ahead(record, &self.ranks, here)?
            {
                self.awaited.entry(mate).or_default().insert(name);
            }
        }
        Ok(Some((contig, templates)))
    }
}

/// Where a read's mate starts, when it is binned at or after the read: `None` for an unpaired
/// read, or a read whose mate is unmapped, is on a contig that is not binned, or comes before it.
fn mate_ahead(
    record: &bam::Record,
    ranks: &HashMap<usize, usize>,
    here: Locus,
) -> PyResult<Option<Locus>> {
    let flags = record.flags();
    if !flags.is_segmented() || flags.is_mate_unmapped() {
        return Ok(None);
    }
    let contig = record
        .mate_reference_sequence_id()
        .transpose()
        .map_err(|error| PyValueError::new_err(error.to_string()))?;
    let start = start_of(record.mate_alignment_start())?;
    let (Some(contig), Some(start)) = (contig, start) else {
        return Ok(None);
    };
    let mate = ranks.get(&contig).map(|&rank| (rank, start));
    Ok(mate.filter(|&mate| mate >= here))
}

/// A 1-based record position as a 0-based one.
fn start_of<P: Into<usize>>(position: Option<io::Result<P>>) -> PyResult<Option<usize>> {
    position
        .transpose()
        .map(|position| position.map(|position| position.into() - 1))
        .map_err(|error| PyValueError::new_err(error.to_string()))
}
