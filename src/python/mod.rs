//! The `streampile._native` extension module: the pileups and tabulation of the Python package,
//! run in Rust over records read by pysam.

// The doc comments of this module are Python docstrings, which name arguments without backticks.
#![allow(clippy::doc_markdown)]

mod bridge;
mod builder;
mod pileup;
mod tabulate;

use std::io;

use pyo3::exceptions::{PyOSError, PyOverflowError, PyValueError};
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

/// A Python `int` given for an option, with its text, kept as a number when it fits an `i64`.
pub(crate) struct Int {
    value: Option<i64>,
    text: String,
}

impl From<i64> for Int {
    fn from(value: i64) -> Self {
        Self {
            value: Some(value),
            text: value.to_string(),
        }
    }
}

impl<'py> FromPyObject<'_, 'py> for Int {
    type Error = PyErr;

    fn extract(object: Borrowed<'_, 'py, PyAny>) -> PyResult<Self> {
        let value = match object.extract::<i64>() {
            Ok(value) => Some(value),
            Err(error) if error.is_instance_of::<PyOverflowError>(object.py()) => None,
            Err(error) => return Err(error),
        };
        Ok(Self {
            value,
            text: object.str()?.to_string(),
        })
    }
}

impl Int {
    /// The option's value as a `T`, or a `ValueError` naming the option when it is not one.
    pub(crate) fn of<T: TryFrom<i64> + Into<i64> + Copy>(self, name: &str, most: T) -> PyResult<T> {
        self.value
            .and_then(|value| T::try_from(value).ok())
            .ok_or_else(|| {
                PyValueError::new_err(format!(
                    "{name} must be from 0 to {}, found: {}",
                    most.into(),
                    self.text
                ))
            })
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

// The bridge reads each pysam record in place, which only the GIL keeps from changing meanwhile.
#[pymodule(gil_used = true)]
fn _native(module: &Bound<'_, PyModule>) -> PyResult<()> {
    bridge::verify(module.py())?;
    module.add_function(wrap_pyfunction!(bridge::direct_bridge, module)?)?;
    module.add_class::<builder::StreamingPileupBuilder>()?;
    module.add_class::<builder::Columns>()?;
    module.add_class::<pileup::Pileup>()?;
    module.add_class::<pileup::PileupRead>()?;
    module.add_class::<pileup::PileupTemplate>()?;
    module.add_class::<tabulate::Tabulation>()?;
    module.add_class::<tabulate::ContigRows>()?;
    module.add_function(wrap_pyfunction!(tabulate::normalize_allele, module)?)?;
    module.add("DEFAULT_EXCLUDE_FLAGS", DEFAULT_EXCLUDE_FLAGS.bits())?;
    module.add("DEFAULT_MIN_BASE_QUALITY", DEFAULT_MIN_BASE_QUALITY)?;
    Ok(())
}
