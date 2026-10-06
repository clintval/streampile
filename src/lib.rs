#![doc = include_str!("../README.crate.md")]

mod auxiliary;
mod builder;
mod error;
mod footprint;
mod overlap;
mod pileup;
#[cfg(feature = "python")]
mod python;
mod source;
mod template;
#[cfg(test)]
mod tests;

pub use auxiliary::{ArraySubtype, AuxArray, AuxElement, AuxValue};
pub use builder::{
    Columns, DEFAULT_EXCLUDE_FLAGS, DEFAULT_MIN_BASE_QUALITY, StreamingPileupBuilder,
};
pub use error::{Error, Result};
pub use overlap::{AgreementStrategy, DisagreementStrategy};
pub use pileup::{EntryKind, MISSING_BASE_QUALITY, Pileup, PileupEntry, PileupTemplate};
pub use source::{AlignmentRecord, RecordSource, Records};
pub use template::{five_prime_distance, template_end_distance};
