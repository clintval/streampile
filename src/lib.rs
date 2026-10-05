//! Forward-only pileups streamed from coordinate-sorted BAM records.
//!
//! A [`StreamingPileupBuilder`] reads records once, from an indexed file or an unindexed pipe,
//! and piles them up at the 0-based positions asked for, moving forward only. Each record's
//! CIGAR is decoded once, when it is read, into a footprint that locates it at every position it
//! spans, and a pileup borrows the records it holds, so building one allocates nothing per entry
//! and per-base auxiliary arrays are read in place.
//!
//! Entries follow the conventions of htslib: a deletion carries the quality of the read's next
//! base, a read with no stored qualities has quality 255 at every base, and an `N` base is a
//! no-call. Insertions are reported at both ends of an alignment.
//!
//! ```no_run
//! use noodles::bam;
//! use streampile::StreamingPileupBuilder;
//!
//! let mut reader = bam::io::reader::Builder::default().build_from_path("reads.bam")?;
//! let header = reader.read_header()?;
//! let mut builder = StreamingPileupBuilder::new(reader, &header)?
//!     .min_mapping_quality(20)
//!     .min_base_quality(30)
//!     .without_overlaps(true);
//!
//! let mut columns = builder.columns("chr1", 0, 1_000)?;
//! while let Some(pileup) = columns.next_pileup() {
//!     let pileup = pileup?;
//!     for entry in pileup.iter().filter(|entry| entry.passes(30)) {
//!         let _ = (entry.base(), entry.is_reverse(), entry.five_prime_distance());
//!     }
//! }
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

mod auxiliary;
mod builder;
mod error;
mod footprint;
mod pileup;
#[cfg(feature = "python")]
mod python;

pub use auxiliary::{ArraySubtype, AuxArray, AuxElement, AuxValue};
pub use builder::{
    Columns, DEFAULT_EXCLUDE_FLAGS, DEFAULT_MIN_BASE_QUALITY, RecordSource, Records,
    StreamingPileupBuilder,
};
pub use error::{Error, Result};
pub use pileup::{EntryKind, MISSING_BASE_QUALITY, Pileup, PileupEntry};
