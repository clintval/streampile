//! The private `streampile._native` extension module, which runs the Rust pileup from Python to
//! compare it with the Python one and to time it.

use std::fs::File;
use std::path::{Path, PathBuf};

use bstr::ByteSlice;
use noodles::bam;
use noodles::bgzf;
use noodles::sam::alignment::record::Flags;
use pyo3::exceptions::{PyOSError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::{EntryKind, Error, Pileup, StreamingPileupBuilder};

type BamReader = bam::io::Reader<bgzf::io::Reader<File>>;

type Entry = (
    String,
    u16,
    &'static str,
    Option<usize>,
    Option<usize>,
    Option<usize>,
    usize,
    Option<String>,
    Option<u8>,
    Option<String>,
    Option<Vec<u8>>,
);

type Column = (Vec<Entry>, usize, usize, String, Vec<u8>);

fn to_python(error: Error) -> PyErr {
    match error {
        Error::Io(error) => PyOSError::new_err(error.to_string()),
        error => PyValueError::new_err(error.to_string()),
    }
}

fn open(path: &Path) -> PyResult<(BamReader, noodles::sam::Header)> {
    let mut reader = bam::io::reader::Builder
        .build_from_path(path)
        .map_err(|error| PyOSError::new_err(format!("{}: {error}", path.display())))?;
    let header = reader
        .read_header()
        .map_err(|error| PyOSError::new_err(error.to_string()))?;
    Ok((reader, header))
}

fn entries(pileup: &Pileup<'_>) -> Vec<Entry> {
    pileup
        .iter()
        .map(|entry| {
            let name = entry
                .record()
                .name()
                .map_or_else(|| "*".to_owned(), ToString::to_string);
            (
                name,
                entry.flags().bits(),
                entry.kind().as_str(),
                entry.query_position(),
                entry.query_position_or_next(),
                entry.insertion_offset(),
                entry.insertion_length(),
                entry.base().map(|base| char::from(base).to_string()),
                entry.quality(),
                entry
                    .inserted_bases()
                    .map(|bases| bases.map(char::from).collect()),
                entry.inserted_qualities().map(Iterator::collect),
            )
        })
        .collect()
}

/// Sweeps spans of a BAM with the Rust builder, in order, and returns every column as its entries,
/// unfiltered depth, filtered depth, bases, and qualities.
#[pyfunction]
#[allow(clippy::needless_pass_by_value, clippy::too_many_arguments)]
#[pyo3(signature = (
    path,
    spans,
    *,
    min_mapping_quality = 0,
    exclude_flags = 0xF00,
    min_base_quality = 13,
    proper_pairs_only = false,
    without_overlaps = false,
))]
fn sweep(
    py: Python<'_>,
    path: PathBuf,
    spans: Vec<(String, usize, usize)>,
    min_mapping_quality: u8,
    exclude_flags: u16,
    min_base_quality: u8,
    proper_pairs_only: bool,
    without_overlaps: bool,
) -> PyResult<Vec<Column>> {
    py.detach(|| {
        let (reader, header) = open(&path)?;
        let mut builder = StreamingPileupBuilder::new(reader, &header)
            .map_err(to_python)?
            .min_mapping_quality(min_mapping_quality)
            .exclude_flags(Flags::from_bits_retain(exclude_flags))
            .min_base_quality(min_base_quality)
            .proper_pairs_only(proper_pairs_only)
            .without_overlaps(without_overlaps);
        let mut swept = Vec::new();
        for (contig, start, end) in spans {
            let mut columns = builder.columns(&contig, start, end).map_err(to_python)?;
            while let Some(pileup) = columns.next_pileup() {
                let pileup = pileup.map_err(to_python)?;
                let bases = pileup.bases().map(char::from).collect();
                swept.push((
                    entries(&pileup),
                    pileup.unfiltered_depth(),
                    pileup.filtered_depth(),
                    bases,
                    pileup.qualities().collect(),
                ));
            }
        }
        builder.close().map_err(to_python)?;
        Ok(swept)
    })
}

/// Counts, at every position of spans of a BAM swept in order, the reads holding each base at a
/// quality floor, the reads with a deletion there (`-`), and the reads with an insertion after it
/// (`+`), as the benchmark's Python engines do.
#[pyclass(module = "streampile._native", unsendable)]
struct ColumnCounts {
    builder: StreamingPileupBuilder<'static, BamReader>,
    spans: Vec<(usize, usize, usize)>,
    span: usize,
    position: usize,
    quality_floor: u8,
}

#[pymethods]
impl ColumnCounts {
    #[new]
    #[allow(clippy::needless_pass_by_value)]
    #[pyo3(signature = (path, spans, *, min_mapping_quality, exclude_flags, quality_floor))]
    fn new(
        path: PathBuf,
        spans: Vec<(String, usize, usize)>,
        min_mapping_quality: u8,
        exclude_flags: u16,
        quality_floor: u8,
    ) -> PyResult<Self> {
        let (reader, header) = open(&path)?;
        let spans = spans
            .into_iter()
            .map(|(contig, start, end)| {
                header
                    .reference_sequences()
                    .get_index_of(contig.as_bytes().as_bstr())
                    .map(|id| (id, start, end))
                    .ok_or_else(|| to_python(Error::UnknownContig(contig)))
            })
            .collect::<PyResult<_>>()?;
        let builder = StreamingPileupBuilder::new(reader, &header)
            .map_err(to_python)?
            .min_mapping_quality(min_mapping_quality)
            .exclude_flags(Flags::from_bits_retain(exclude_flags));
        let position = 0;
        Ok(Self {
            builder,
            spans,
            span: 0,
            position,
            quality_floor,
        })
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__<'py>(&mut self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyDict>>> {
        let (reference_sequence_id, position) = loop {
            let Some(&(id, start, end)) = self.spans.get(self.span) else {
                return Ok(None);
            };
            let position = self.position.max(start);
            if position < end {
                self.position = position + 1;
                break (id, position);
            }
            self.span += 1;
            self.position = 0;
        };
        let pileup = self
            .builder
            .pileup_at(reference_sequence_id, position)
            .map_err(to_python)?;
        let mut counts = [0_u32; 256];
        for entry in pileup.iter() {
            if entry.kind() == EntryKind::Insertion {
                counts[usize::from(b'+')] += 1;
            } else if entry.passes(self.quality_floor) {
                let key = if entry.is_deletion() {
                    b'-'
                } else {
                    entry.base().unwrap_or(b'N')
                };
                counts[usize::from(key)] += 1;
            }
        }
        let column = PyDict::new(py);
        for (key, &count) in counts.iter().enumerate() {
            if count > 0 {
                column.set_item(char::from(key as u8), count)?;
            }
        }
        Ok(Some(column))
    }
}

#[pymodule(gil_used = false)]
fn _native(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(sweep, module)?)?;
    module.add_class::<ColumnCounts>()?;
    Ok(())
}
