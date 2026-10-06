//! Copies a pysam `AlignedSegment` into a BAM record.
//!
//! An `AlignedSegment` holds its record as htslib's `bam1_t`, behind the `_delegate` pointer that
//! pysam's `libcalignedsegment.pxd` declares first in the object, after the Python object header
//! and Cython's vtable pointer. Reading the record from there takes one copy, where reading it
//! through the object's attributes takes a dozen Python calls. The layout is checked once, against
//! a record made from known SAM text, and if it differs, records are read through attributes.

use std::cell::RefCell;
use std::ffi::c_void;
use std::mem::size_of;
use std::sync::atomic::{AtomicBool, Ordering};

use noodles::bam;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::sync::PyOnceLock;
use pyo3::types::{PyString, PyType};

const DELEGATE_OFFSET: usize = size_of::<pyo3::ffi::PyObject>() + size_of::<*const c_void>();

const MOST_CIGAR_OPERATORS: usize = 0xFFFF;

const SEQUENCE_CODES: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

static VERIFIED: AtomicBool = AtomicBool::new(false);
static DIRECT: AtomicBool = AtomicBool::new(false);

static ALIGNED_SEGMENT: PyOnceLock<Py<PyType>> = PyOnceLock::new();

thread_local! {
    static STAGED: RefCell<Vec<u8>> = const { RefCell::new(Vec::new()) };
}

/// htslib's `bam1_core_t`, as htslib 1.10 and later lay it out.
#[repr(C)]
struct Core {
    pos: i64,
    tid: i32,
    bin: u16,
    qual: u8,
    l_extranul: u8,
    flag: u16,
    l_qname: u16,
    n_cigar: u32,
    l_qseq: i32,
    mtid: i32,
    mpos: i64,
    isize: i64,
}

/// htslib's `bam1_t`, up to the fields read here.
#[repr(C)]
struct Delegate {
    core: Core,
    id: u64,
    data: *const u8,
    l_data: i32,
}

/// The fixed-length fields of a BAM record, before its variable-length ones.
struct Head {
    tid: i32,
    pos: i64,
    mapq: u8,
    bin: u16,
    flag: u16,
    l_seq: usize,
    mtid: i32,
    mpos: i64,
    tlen: i64,
}

/// Reads a pysam `AlignedSegment`, or any object with its attributes, into a BAM record.
///
/// The object must not change while this runs, which holding the GIL and a reference to it
/// ensures; nothing of it is kept, written, or freed.
pub(crate) fn read(object: &Bound<'_, PyAny>, record: &mut bam::Record) -> PyResult<()> {
    STAGED.with_borrow_mut(|staged| {
        staged.clear();
        if DIRECT.load(Ordering::Relaxed) && is_aligned_segment(object)? {
            // SAFETY: the object is an `AlignedSegment`, whose layout was checked at import.
            unsafe { stage_delegate(object, staged) }?;
        } else {
            stage_attributes(object, staged)?;
        }
        bam::io::Reader::from(&staged[..])
            .read_record(record)
            .map_err(|error| PyValueError::new_err(format!("invalid record: {error}")))?;
        Ok(())
    })
}

/// Whether records are read from htslib's `bam1_t`, and with `enabled`, sets it first: records
/// are read from `bam1_t` only when its layout was found as expected.
#[pyfunction]
#[pyo3(signature = (enabled = None))]
pub(crate) fn direct_bridge(enabled: Option<bool>) -> bool {
    if let Some(enabled) = enabled {
        DIRECT.store(
            enabled && VERIFIED.load(Ordering::Relaxed),
            Ordering::Relaxed,
        );
    }
    DIRECT.load(Ordering::Relaxed)
}

/// Checks pysam's layout against a record made from known SAM text, read through its
/// attributes, and reads records from `bam1_t` from now on only if both agree.
pub(crate) fn verify(py: Python<'_>) -> PyResult<()> {
    const HEADER: &str = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n";
    const LINE: &str = "pair/1\t99\tchr1\t11\t42\t3S4M1I2M2D1M\t=\t20\t30\tACGTACGTACG\t\
        ABCDEFGHIJK\tMC:Z:10M\tXI:i:-7\tXB:B:s,1,-2";
    let pysam = py.import("pysam")?;
    let header = pysam
        .getattr("AlignmentHeader")?
        .call_method1("from_text", (HEADER,))?;
    let segment = pysam
        .getattr("AlignedSegment")?
        .call_method1("fromstring", (LINE, header))?;
    let basic_size: usize = segment.get_type().getattr("__basicsize__")?.extract()?;
    let mut through_attributes = bam::Record::default();
    STAGED.with_borrow_mut(|staged| -> PyResult<()> {
        staged.clear();
        stage_attributes(&segment, staged)?;
        bam::io::Reader::from(&staged[..])
            .read_record(&mut through_attributes)
            .map_err(|error| PyValueError::new_err(format!("invalid record: {error}")))?;
        Ok(())
    })?;
    let mut direct = bam::Record::default();
    let agrees = basic_size >= DELEGATE_OFFSET + size_of::<*const c_void>()
        && STAGED.with_borrow_mut(|staged| {
            staged.clear();
            // SAFETY: the object is large enough to hold the pointer read, and the pointer read is
            // checked for null; what it points to is compared with the attributes below.
            unsafe { stage_delegate(&segment, staged) }.is_ok()
                && bam::io::Reader::from(&staged[..])
                    .read_record(&mut direct)
                    .is_ok()
        })
        && same_record(&direct, &through_attributes)
        && direct.data().as_bytes() == b"MCZ10M\0XIc\xf9XBBs\x02\0\0\0\x01\0\xfe\xff";
    VERIFIED.store(agrees, Ordering::Relaxed);
    DIRECT.store(agrees, Ordering::Relaxed);
    Ok(())
}

fn same_record(a: &bam::Record, b: &bam::Record) -> bool {
    a.reference_sequence_id().transpose().ok() == b.reference_sequence_id().transpose().ok()
        && a.alignment_start().transpose().ok() == b.alignment_start().transpose().ok()
        && a.mapping_quality() == b.mapping_quality()
        && a.flags() == b.flags()
        && a.mate_reference_sequence_id().transpose().ok()
            == b.mate_reference_sequence_id().transpose().ok()
        && a.mate_alignment_start().transpose().ok() == b.mate_alignment_start().transpose().ok()
        && a.template_length() == b.template_length()
        && a.name() == b.name()
        && a.cigar().as_bytes() == b.cigar().as_bytes()
        && a.sequence().as_bytes() == b.sequence().as_bytes()
        && a.quality_scores().as_bytes() == b.quality_scores().as_bytes()
        && a.name().is_some_and(|name| name == "pair/1")
}

fn is_aligned_segment(object: &Bound<'_, PyAny>) -> PyResult<bool> {
    let py = object.py();
    let aligned_segment = ALIGNED_SEGMENT.get_or_try_init(py, || {
        py.import("pysam")?
            .getattr("AlignedSegment")?
            .cast_into::<PyType>()
            .map(Bound::unbind)
            .map_err(PyErr::from)
    })?;
    let kind = object.get_type();
    Ok(kind.is(aligned_segment.bind(py)) || kind.is_subclass(aligned_segment.bind(py))?)
}

/// Stages the record of an `AlignedSegment` from its `bam1_t`.
///
/// # Safety
///
/// `object` must be a pysam `AlignedSegment` laid out as [`verify`] checks.
unsafe fn stage_delegate(object: &Bound<'_, PyAny>, staged: &mut Vec<u8>) -> PyResult<()> {
    let invalid = || PyValueError::new_err("an AlignedSegment holds no valid record");
    // SAFETY: the caller guarantees `_delegate` is at this offset of the object.
    let delegate = unsafe {
        object
            .as_ptr()
            .cast::<*const Delegate>()
            .byte_add(DELEGATE_OFFSET)
            .read()
    };
    // SAFETY: a non-null `_delegate` points at the `bam1_t` the object owns.
    let delegate = unsafe { delegate.as_ref() }.ok_or_else(invalid)?;
    let core = &delegate.core;
    let length = usize::try_from(delegate.l_data).map_err(|_| invalid())?;
    let data = if length == 0 || delegate.data.is_null() {
        &[][..]
    } else {
        // SAFETY: htslib keeps `l_data` bytes of the record at `data`.
        unsafe { std::slice::from_raw_parts(delegate.data, length) }
    };
    let l_qname = usize::from(core.l_qname);
    let l_seq = usize::try_from(core.l_qseq).map_err(|_| invalid())?;
    let operators = core.n_cigar as usize;
    let name_end = l_qname
        .checked_sub(usize::from(core.l_extranul))
        .ok_or_else(invalid)?;
    let cigar_end = l_qname + operators * 4;
    let sequence_end = cigar_end + l_seq.div_ceil(2);
    let qualities_end = sequence_end + l_seq;
    if qualities_end > data.len() {
        return Err(invalid());
    }
    let name = if name_end == 0 {
        b"*\0"
    } else {
        &data[..name_end]
    };
    let head = Head {
        tid: core.tid,
        pos: core.pos,
        mapq: core.qual,
        bin: core.bin,
        flag: core.flag,
        l_seq,
        mtid: core.mtid,
        mpos: core.mpos,
        tlen: core.isize,
    };
    let cigar = &data[l_qname..cigar_end];
    let ops = cigar
        .as_chunks::<4>()
        .0
        .iter()
        .map(|op| u32::from_ne_bytes(*op));
    stage(
        staged,
        &head,
        name,
        ops,
        operators,
        |staged| staged.extend_from_slice(&data[cigar_end..qualities_end]),
        &data[qualities_end..],
    )
}

/// Stages the record of any object with the attributes of an `AlignedSegment`, keeping, of its
/// tags, only the `MC` tag a pileup reads.
fn stage_attributes(object: &Bound<'_, PyAny>, staged: &mut Vec<u8>) -> PyResult<()> {
    let name: Option<String> = object.getattr("query_name")?.extract()?;
    let cigar: Option<Vec<(u32, u32)>> = object.getattr("cigartuples")?.extract()?;
    let sequence: Option<String> = object.getattr("query_sequence")?.extract()?;
    let qualities: Option<Vec<u8>> = object.getattr("query_qualities")?.extract()?;
    let sequence = sequence.unwrap_or_default();
    let head = Head {
        tid: object.getattr("reference_id")?.extract()?,
        pos: object.getattr("reference_start")?.extract()?,
        mapq: object.getattr("mapping_quality")?.extract()?,
        bin: 4680,
        flag: object.getattr("flag")?.extract()?,
        l_seq: sequence.len(),
        mtid: object.getattr("next_reference_id")?.extract()?,
        mpos: object.getattr("next_reference_start")?.extract()?,
        tlen: object.getattr("template_length")?.extract()?,
    };
    let mut aux = Vec::new();
    if object.call_method1("has_tag", ("MC",))?.is_truthy()? {
        let (value, kind): (Bound<'_, PyAny>, String) =
            object.call_method1("get_tag", ("MC", true))?.extract()?;
        aux.extend_from_slice(b"MC");
        stage_value(&value, kind.as_bytes(), &mut aux)?;
    }
    let mut name = name.unwrap_or_else(|| "*".to_owned()).into_bytes();
    name.push(0);
    let cigar = cigar.unwrap_or_default();
    let operators = cigar.len();
    let ops = cigar.into_iter().map(|(kind, length)| length << 4 | kind);
    let quality_length = head.l_seq;
    stage(
        staged,
        &head,
        &name,
        ops,
        operators,
        |staged| {
            let bases = sequence.as_bytes();
            for pair in bases.chunks(2) {
                let high = nibble(pair[0]) << 4;
                staged.push(high | pair.get(1).map_or(0, |&base| nibble(base)));
            }
            match qualities.as_deref() {
                Some(qualities) if qualities.len() == quality_length => {
                    staged.extend_from_slice(qualities);
                }
                _ => staged.extend(std::iter::repeat_n(0xFF, quality_length)),
            }
        },
        &aux,
    )
}

/// Stages a tag's type and value as BAM stores them, from the value and the type pysam reads.
fn stage_value(value: &Bound<'_, PyAny>, kind: &[u8], staged: &mut Vec<u8>) -> PyResult<()> {
    match kind {
        [b'Z' | b'H'] => {
            staged.extend_from_slice(kind);
            staged.extend_from_slice(value.cast::<PyString>()?.to_str()?.as_bytes());
            staged.push(0);
        }
        [b'A'] => {
            staged.extend_from_slice(kind);
            let text = value.cast::<PyString>()?.to_str()?;
            staged.push(text.bytes().next().unwrap_or(b' '));
        }
        [b'B', subtype] => {
            staged.extend_from_slice(kind);
            let length = u32::try_from(value.len()?)
                .map_err(|_| PyValueError::new_err("An array tag is too long."))?;
            staged.extend_from_slice(&length.to_le_bytes());
            for element in value.try_iter()? {
                stage_number(&element?, *subtype, staged)?;
            }
        }
        [number] => {
            staged.push(*number);
            stage_number(value, *number, staged)?;
        }
        _ => {
            return Err(PyValueError::new_err(format!(
                "A tag has the unknown type {}.",
                String::from_utf8_lossy(kind)
            )));
        }
    }
    Ok(())
}

/// Stages a number of one of BAM's numeric types, little-endian.
fn stage_number(value: &Bound<'_, PyAny>, kind: u8, staged: &mut Vec<u8>) -> PyResult<()> {
    match kind {
        b'c' => staged.extend_from_slice(&value.extract::<i8>()?.to_le_bytes()),
        b'C' => staged.extend_from_slice(&value.extract::<u8>()?.to_le_bytes()),
        b's' => staged.extend_from_slice(&value.extract::<i16>()?.to_le_bytes()),
        b'S' => staged.extend_from_slice(&value.extract::<u16>()?.to_le_bytes()),
        b'i' => staged.extend_from_slice(&value.extract::<i32>()?.to_le_bytes()),
        b'I' => staged.extend_from_slice(&value.extract::<u32>()?.to_le_bytes()),
        b'f' => staged.extend_from_slice(&value.extract::<f32>()?.to_le_bytes()),
        b'd' => staged.extend_from_slice(&value.extract::<f64>()?.to_le_bytes()),
        other => {
            return Err(PyValueError::new_err(format!(
                "A tag has the unknown type {}.",
                char::from(other)
            )));
        }
    }
    Ok(())
}

fn nibble(base: u8) -> u8 {
    let base = base.to_ascii_uppercase();
    SEQUENCE_CODES
        .iter()
        .position(|&code| code == base)
        .map_or(15, |code| code as u8)
}

/// Stages a record as a BAM reader reads it, block size first, moving a CIGAR of more operators
/// than BAM holds to a `CG` tag, as htslib writes it.
fn stage(
    staged: &mut Vec<u8>,
    head: &Head,
    name: &[u8],
    ops: impl Iterator<Item = u32> + Clone,
    operators: usize,
    sequence_and_qualities: impl FnOnce(&mut Vec<u8>),
    aux: &[u8],
) -> PyResult<()> {
    let too_long = |field: &str| PyValueError::new_err(format!("a record's {field} is too long"));
    let narrow = |value: i64, field: &str| i32::try_from(value).map_err(|_| too_long(field));
    let l_read_name = u8::try_from(name.len()).map_err(|_| too_long("name"))?;
    let long_cigar = operators > MOST_CIGAR_OPERATORS;
    staged.extend_from_slice(&[0; 4]);
    staged.extend_from_slice(&head.tid.to_le_bytes());
    staged.extend_from_slice(&narrow(head.pos, "position")?.to_le_bytes());
    staged.push(l_read_name);
    staged.push(head.mapq);
    staged.extend_from_slice(&head.bin.to_le_bytes());
    let stored_operators = if long_cigar { 2 } else { operators as u16 };
    staged.extend_from_slice(&stored_operators.to_le_bytes());
    staged.extend_from_slice(&head.flag.to_le_bytes());
    let l_seq = u32::try_from(head.l_seq).map_err(|_| too_long("sequence"))?;
    staged.extend_from_slice(&l_seq.to_le_bytes());
    staged.extend_from_slice(&head.mtid.to_le_bytes());
    staged.extend_from_slice(&narrow(head.mpos, "mate position")?.to_le_bytes());
    staged.extend_from_slice(&narrow(head.tlen, "template length")?.to_le_bytes());
    staged.extend_from_slice(name);
    if long_cigar {
        let span: u32 = ops
            .clone()
            .filter(|op| matches!(op & 0xF, 0 | 2 | 3 | 7 | 8))
            .map(|op| op >> 4)
            .sum();
        staged.extend_from_slice(&(l_seq << 4 | 4).to_le_bytes());
        staged.extend_from_slice(&(span << 4 | 3).to_le_bytes());
    } else {
        for op in ops.clone() {
            staged.extend_from_slice(&op.to_le_bytes());
        }
    }
    sequence_and_qualities(staged);
    staged.extend_from_slice(aux);
    if long_cigar {
        staged.extend_from_slice(b"CGBI");
        staged.extend_from_slice(&(operators as u32).to_le_bytes());
        for op in ops {
            staged.extend_from_slice(&op.to_le_bytes());
        }
    }
    let block_size = u32::try_from(staged.len() - 4).map_err(|_| too_long("record"))?;
    staged[..4].copy_from_slice(&block_size.to_le_bytes());
    Ok(())
}
