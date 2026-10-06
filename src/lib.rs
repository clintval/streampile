#![doc = include_str!("../README.crate.md")]

mod auxiliary;
mod builder;
mod error;
mod footprint;
mod pileup;
#[cfg(feature = "python")]
mod python;
mod source;
#[cfg(test)]
mod tests;

pub use auxiliary::{ArraySubtype, AuxArray, AuxElement, AuxValue};
pub use builder::{
    Columns, DEFAULT_EXCLUDE_FLAGS, DEFAULT_MIN_BASE_QUALITY, StreamingPileupBuilder,
};
pub use error::{Error, Result};
pub use pileup::{EntryKind, MISSING_BASE_QUALITY, Pileup, PileupEntry};
pub use source::{AlignmentRecord, RecordSource, Records};
