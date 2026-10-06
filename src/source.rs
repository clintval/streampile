use std::io::{self, Read};

use noodles::bam;

/// A record a builder piles up: a BAM record, or a BAM record with something its source keeps
/// beside it, such as the object it was read from.
pub trait AlignmentRecord: Default {
    /// The BAM record.
    fn bam(&self) -> &bam::Record;

    /// Lets go of anything kept beside the BAM record once a builder is done with the record.
    ///
    /// A builder fills a released record again with a later one, so a BAM record keeps its
    /// buffers, and by default nothing is let go of.
    fn release(&mut self) {}
}

impl AlignmentRecord for bam::Record {
    fn bam(&self) -> &bam::Record {
        self
    }
}

/// A coordinate-sorted stream of records, which may be an unindexed pipe.
pub trait RecordSource {
    /// The records the source reads: [`bam::Record`] for a BAM reader.
    type Record: AlignmentRecord;

    /// Reads the next record into `record`, reusing its buffers, and returns `false` at the end.
    fn read_record(&mut self, record: &mut Self::Record) -> io::Result<bool>;
}

impl<R: Read> RecordSource for bam::io::Reader<R> {
    type Record = bam::Record;

    fn read_record(&mut self, record: &mut bam::Record) -> io::Result<bool> {
        bam::io::Reader::read_record(self, record).map(|read| read > 0)
    }
}

impl<S: RecordSource + ?Sized> RecordSource for &mut S {
    type Record = S::Record;

    fn read_record(&mut self, record: &mut S::Record) -> io::Result<bool> {
        (**self).read_record(record)
    }
}

impl<S: RecordSource + ?Sized> RecordSource for Box<S> {
    type Record = S::Record;

    fn read_record(&mut self, record: &mut S::Record) -> io::Result<bool> {
        (**self).read_record(record)
    }
}

/// A [`RecordSource`] over any iterator of records.
#[derive(Debug)]
pub struct Records<I>(I);

impl<I> Records<I> {
    /// Wraps an iterator of records.
    pub fn new(records: I) -> Self {
        Self(records)
    }
}

impl<R: AlignmentRecord, I: Iterator<Item = io::Result<R>>> RecordSource for Records<I> {
    type Record = R;

    fn read_record(&mut self, record: &mut R) -> io::Result<bool> {
        match self.0.next() {
            Some(next) => {
                *record = next?;
                Ok(true)
            }
            None => Ok(false),
        }
    }
}
