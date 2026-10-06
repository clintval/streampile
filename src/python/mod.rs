//! The `streampile._native` extension module: the pileups and tabulation of the Python package,
//! run in Rust over records read by pysam.

// The doc comments of this module are Python docstrings, which name arguments without backticks.
#![allow(clippy::doc_markdown)]

mod bridge;
mod builder;
mod pileup;
mod tabulate;

use std::io;

use pyo3::exceptions::{PyOSError, PyValueError};
use pyo3::prelude::*;

use crate::{DEFAULT_EXCLUDE_FLAGS, DEFAULT_MIN_BASE_QUALITY, Error};

/// The Python exception of an error: the very exception a Python callback raised, or a
/// `ValueError` worded as the Python package words it.
fn to_python(error: Error) -> PyErr {
    match error {
        Error::Io(error) => from_io(error),
        Error::NotCoordinateSorted { .. } => {
            PyValueError::new_err("Records must be coordinate sorted.")
        }
        Error::OutOfOrder { name } => {
            PyValueError::new_err(format!("Records are out of coordinate order at {name}."))
        }
        Error::UnknownContig(contig) => {
            PyValueError::new_err(format!("Contig {contig} is not in the header."))
        }
        Error::Closed => PyValueError::new_err("The builder is closed."),
        error => PyValueError::new_err(sentence(&error.to_string())),
    }
}

fn from_io(error: io::Error) -> PyErr {
    if let Some(inner) = error.get_ref()
        && inner.is::<PyErr>()
    {
        let inner = error.into_inner().expect("an error inside");
        return *inner.downcast::<PyErr>().expect("a Python exception");
    }
    match error.kind() {
        io::ErrorKind::InvalidData | io::ErrorKind::UnexpectedEof => {
            PyValueError::new_err(sentence(&error.to_string()))
        }
        _ => PyOSError::new_err(error.to_string()),
    }
}

/// A message as a sentence: capitalized, and ending with a period.
fn sentence(message: &str) -> String {
    let mut characters = message.chars();
    let mut sentence: String = characters
        .next()
        .map(|first| first.to_uppercase().chain(characters).collect())
        .unwrap_or_default();
    if !sentence.ends_with(['.', '!', '?']) {
        sentence.push('.');
    }
    sentence
}

#[pymodule(gil_used = false)]
fn _native(module: &Bound<'_, PyModule>) -> PyResult<()> {
    bridge::verify(module.py())?;
    module.add_function(wrap_pyfunction!(bridge::direct_bridge, module)?)?;
    module.add_class::<builder::StreamingPileupBuilder>()?;
    module.add_class::<pileup::Pileup>()?;
    module.add_class::<pileup::PileupRead>()?;
    module.add_class::<pileup::PileupTemplate>()?;
    module.add_class::<tabulate::Tabulation>()?;
    module.add_function(wrap_pyfunction!(tabulate::normalize_allele, module)?)?;
    module.add("DEFAULT_EXCLUDE_FLAGS", DEFAULT_EXCLUDE_FLAGS.bits())?;
    module.add("DEFAULT_MIN_BASE_QUALITY", DEFAULT_MIN_BASE_QUALITY)?;
    Ok(())
}
