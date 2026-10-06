use std::io;

/// An error from building pileups.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum Error {
    /// Reading a record failed, or a record is malformed.
    #[error(transparent)]
    Io(#[from] io::Error),

    /// The header does not declare `SO:coordinate`.
    #[error("records must be coordinate sorted, but the header declares {}", found.as_deref().unwrap_or("no sort order"))]
    NotCoordinateSorted {
        /// The sort order the header declares, if any.
        found: Option<String>,
    },

    /// A record is malformed, such as one whose CIGAR and sequence differ in length.
    #[error("read {name} is invalid: {source}")]
    InvalidRecord {
        /// The name of the record.
        name: String,
        /// What is wrong with it.
        source: io::Error,
    },

    /// A forward read of an FR pair has no `MC` tag to find its mate's 5′ end with.
    #[error("read {name} has no MC tag to find its mate's 5' end with")]
    MissingMateCigar {
        /// The name of the read.
        name: String,
    },

    /// A forward read of an FR pair has an `MC` tag that is not a CIGAR string spanning at least
    /// one reference base.
    #[error("read {name} has an invalid MC tag: {value}")]
    InvalidMateCigar {
        /// The name of the read.
        name: String,
        /// The value of the tag.
        value: String,
    },

    /// A record starts before the record read just before it.
    #[error("records are out of coordinate order at {name}")]
    OutOfOrder {
        /// The name of the first record out of order.
        name: String,
    },

    /// A pileup was asked for before the last one asked for.
    #[error("attempted to advance to {contig}:{position} from {from_contig}:{from_position}")]
    Backwards {
        /// The contig asked for.
        contig: String,
        /// The 0-based position asked for.
        position: usize,
        /// The contig of the last pileup.
        from_contig: String,
        /// The 0-based position of the last pileup.
        from_position: usize,
    },

    /// A contig name is not in the header.
    #[error("contig {0} is not in the header")]
    UnknownContig(String),

    /// A reference sequence ID is not in the header.
    #[error("reference sequence ID {0} is not in the header")]
    UnknownReferenceSequenceId(usize),

    /// A span ends before it starts.
    #[error("end {end} is before start {start}")]
    InvalidSpan {
        /// The 0-based first position.
        start: usize,
        /// The 0-based position after the last one.
        end: usize,
    },

    /// The builder was asked for a pileup after it was closed.
    #[error("the builder is closed")]
    Closed,

    /// A name is not the name of any strategy of a kind.
    #[error("'{name}' is not a valid {kind}")]
    UnknownStrategy {
        /// The kind of strategy, such as `AgreementStrategy`.
        kind: &'static str,
        /// The name.
        name: String,
    },

    /// The builder was asked for a pileup after an error stopped it while advancing.
    #[error("the builder stopped at an earlier error")]
    Stopped,
}

/// A result whose error is an [`Error`].
pub type Result<T> = std::result::Result<T, Error>;

/// An error of malformed data, as noodles reports one.
pub(crate) fn invalid_data(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.into())
}
