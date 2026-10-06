mod builder;
mod fgbio;
mod fixture;
mod footprint;
mod orientation;
mod pileup;
mod templates;
mod testing;

use std::io;

use noodles::bam;
use noodles::sam::{self, alignment::io::Write as _};

use crate::{Pileup, Records, StreamingPileupBuilder};

pub(crate) const HEADER: &str =
    "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:1000\n";

/// A read built field by field, with Q40 bases unless given, and `*` as in SAM.
#[derive(Clone, Debug)]
pub(crate) struct Read {
    name: String,
    flag: u16,
    contig: String,
    start: Option<usize>,
    mapq: u8,
    cigar: String,
    mate: Option<(String, usize)>,
    tlen: i32,
    bases: String,
    quals: Option<Vec<u8>>,
    tags: Vec<String>,
}

/// A mapped read at a 0-based start.
pub(crate) fn read(name: &str, start: usize, cigar: &str, bases: &str) -> Read {
    let quals = (bases != "*").then(|| vec![40; bases.len()]);
    Read {
        name: name.to_owned(),
        flag: 0,
        contig: "chr1".to_owned(),
        start: Some(start),
        mapq: 60,
        cigar: cigar.to_owned(),
        mate: None,
        tlen: 0,
        bases: bases.to_owned(),
        quals,
        tags: Vec::new(),
    }
}

/// A read with no position at all.
pub(crate) fn unmapped(name: &str) -> Read {
    Read {
        flag: 4,
        contig: "*".to_owned(),
        start: None,
        mapq: 0,
        ..read(name, 0, "*", "ACGT")
    }
}

impl Read {
    pub(crate) fn flag(mut self, flag: u16) -> Self {
        self.flag = flag;
        self
    }

    pub(crate) fn mapq(mut self, mapq: u8) -> Self {
        self.mapq = mapq;
        self
    }

    pub(crate) fn contig(mut self, contig: &str) -> Self {
        contig.clone_into(&mut self.contig);
        self
    }

    pub(crate) fn quals(mut self, quals: &[u8]) -> Self {
        self.quals = Some(quals.to_vec());
        self
    }

    pub(crate) fn no_quals(mut self) -> Self {
        self.quals = None;
        self
    }

    pub(crate) fn mate(mut self, contig: &str, start: usize, tlen: i32) -> Self {
        self.mate = Some((contig.to_owned(), start));
        self.tlen = tlen;
        self
    }

    pub(crate) fn tag(mut self, tag: &str) -> Self {
        self.tags.push(tag.to_owned());
        self
    }

    fn to_sam(&self) -> String {
        let position = self.start.map_or(0, |start| start + 1);
        let (rnext, pnext) = match &self.mate {
            Some((contig, start)) if *contig == self.contig => ("=".to_owned(), start + 1),
            Some((contig, start)) => (contig.clone(), start + 1),
            None => ("*".to_owned(), 0),
        };
        let quals = self.quals.as_ref().map_or_else(
            || "*".to_owned(),
            |quals| {
                quals
                    .iter()
                    .map(|&quality| char::from(quality + 33))
                    .collect()
            },
        );
        let mut line = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.name,
            self.flag,
            self.contig,
            position,
            self.mapq,
            self.cigar,
            rnext,
            pnext,
            self.tlen,
            self.bases,
            quals,
        );
        for tag in &self.tags {
            line.push('\t');
            line.push_str(tag);
        }
        line.push('\n');
        line
    }
}

/// Writes reads to an in-memory BAM, BGZF-compressed, under a header.
pub(crate) fn bam_bytes(header: &str, reads: &[Read]) -> Vec<u8> {
    let text: String = std::iter::once(header.to_owned())
        .chain(reads.iter().map(Read::to_sam))
        .collect();
    let mut sam_reader = sam::io::Reader::new(text.as_bytes());
    let header = sam_reader.read_header().expect("a valid header");
    let mut writer = bam::io::Writer::new(Vec::new());
    writer.write_header(&header).expect("the header is written");
    for record in sam_reader.record_bufs(&header) {
        let record = record.expect("a valid SAM record");
        writer
            .write_alignment_record(&header, &record)
            .expect("the record is written");
    }
    writer.into_inner().finish().expect("the BAM is finished")
}

/// The header and records of reads, read back from a BAM.
pub(crate) fn records(header: &str, reads: &[Read]) -> (sam::Header, Vec<bam::Record>) {
    let bytes = bam_bytes(header, reads);
    let mut reader = bam::io::Reader::new(&bytes[..]);
    let header = reader.read_header().expect("a valid header");
    let records = reader
        .records()
        .collect::<io::Result<_>>()
        .expect("valid records");
    (header, records)
}

/// A record on chr1 built from raw BAM fields, which may disagree in ways a BAM writer refuses.
pub(crate) fn raw_record(name: &str, start: i32, cigar: &[u32], bases: usize) -> bam::Record {
    let mut block = Vec::new();
    block.extend_from_slice(&0_i32.to_le_bytes());
    block.extend_from_slice(&start.to_le_bytes());
    block.extend_from_slice(&[(name.len() + 1) as u8, 60]);
    block.extend_from_slice(&0_u16.to_le_bytes());
    block.extend_from_slice(&(cigar.len() as u16).to_le_bytes());
    block.extend_from_slice(&0_u16.to_le_bytes());
    block.extend_from_slice(&(bases as u32).to_le_bytes());
    block.extend_from_slice(&(-1_i32).to_le_bytes());
    block.extend_from_slice(&(-1_i32).to_le_bytes());
    block.extend_from_slice(&0_i32.to_le_bytes());
    block.extend_from_slice(name.as_bytes());
    block.push(0);
    for op in cigar {
        block.extend_from_slice(&op.to_le_bytes());
    }
    block.extend(std::iter::repeat_n(0x12, bases.div_ceil(2)));
    block.extend(std::iter::repeat_n(30, bases));
    let mut bytes = (block.len() as u32).to_le_bytes().to_vec();
    bytes.extend_from_slice(&block);
    let mut record = bam::Record::default();
    bam::io::Reader::from(&bytes[..])
        .read_record(&mut record)
        .expect("a record-like block");
    record
}

pub(crate) type Source = Records<
    std::iter::Map<std::vec::IntoIter<bam::Record>, fn(bam::Record) -> io::Result<bam::Record>>,
>;

/// A source of records that are already read.
pub(crate) fn source(records: Vec<bam::Record>) -> Source {
    Records::new(
        records
            .into_iter()
            .map(Ok as fn(bam::Record) -> io::Result<bam::Record>),
    )
}

/// A builder over reads, already in coordinate order, under the default header.
pub(crate) fn builder<'f>(reads: &[Read]) -> StreamingPileupBuilder<'f, Source> {
    let (header, records) = records(HEADER, reads);
    StreamingPileupBuilder::new(source(records), &header).expect("a coordinate-sorted header")
}

/// A builder that piles up every placed read, as Python's `Pileup.from_alignments` does.
pub(crate) fn unfiltered<'f>(reads: &[Read]) -> StreamingPileupBuilder<'f, Source> {
    builder(reads).exclude_flags(sam::alignment::record::Flags::empty())
}

pub(crate) type Entry = (
    String,
    &'static str,
    Option<usize>,
    Option<usize>,
    Option<String>,
);

/// Each entry of a pileup as (read name, kind, query position, next, inserted bases).
pub(crate) fn entries(pileup: &Pileup<'_>) -> Vec<Entry> {
    pileup
        .iter()
        .map(|entry| {
            (
                name(entry.record()),
                entry.kind().as_str(),
                entry.query_position(),
                entry.query_position_or_next(),
                entry
                    .inserted_bases()
                    .map(|bases| bases.map(char::from).collect()),
            )
        })
        .collect()
}

pub(crate) fn entry(
    name: &str,
    kind: &'static str,
    position: Option<usize>,
    next: Option<usize>,
    inserted: Option<&str>,
) -> Entry {
    (
        name.to_owned(),
        kind,
        position,
        next,
        inserted.map(ToOwned::to_owned),
    )
}

pub(crate) fn name(record: &bam::Record) -> String {
    record
        .name()
        .map_or_else(|| "*".to_owned(), ToString::to_string)
}

pub(crate) fn names(pileup: &Pileup<'_>) -> Vec<String> {
    pileup.iter().map(|entry| name(entry.record())).collect()
}

pub(crate) fn bases(pileup: &Pileup<'_>) -> String {
    pileup.bases().map(char::from).collect()
}
