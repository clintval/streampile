//! Test records built as fgbio's `SamBuilder` builds them, with full mate information on pairs,
//! and a way to pile them up in memory with a [`StreamingPileupBuilder`].
//!
//! Builder fields carry fgbio's parameter names and defaults, and starts are 1-based, with 0 for a
//! read with no position.

use std::fs::File;
use std::io;
use std::num::NonZeroUsize;
use std::path::Path;

use noodles::bam;
use noodles::core::Position;
use noodles::sam;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::io::Write as _;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record::{Flags, MappingQuality};
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::record_buf::{Cigar, Data, QualityScores, Sequence};
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::header::tag::SORT_ORDER;
use noodles::sam::header::record::value::map::read_group::tag::SAMPLE;
use noodles::sam::header::record::value::map::{
    Header, ReadGroup, ReferenceSequence, header::Version,
};

use crate::{Records, StreamingPileupBuilder};

/// The read group of every built record, as in fgbio.
pub const READ_GROUP_ID: &str = "A";

/// The sample of the built records' read group, as in fgbio.
pub const SAMPLE_NAME: &str = "Sample";

/// A source of built records, in coordinate order, that a [`StreamingPileupBuilder`] reads.
pub type BuiltRecords = Records<std::vec::IntoIter<io::Result<bam::Record>>>;

/// The strand a built read aligns to.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Strand {
    /// The forward strand.
    #[default]
    Plus,
    /// The reverse strand.
    Minus,
}

/// A read pair to build, with the defaults of fgbio's `SamBuilder.addPair`: an FR pair, the
/// first read forward.
#[derive(Clone, Debug)]
pub struct Pair {
    /// The name of both reads, or the next sequential name.
    pub name: Option<String>,
    /// The first read's bases, or random bases of the builder's read length.
    pub bases1: Option<String>,
    /// The second read's bases, or random bases of the builder's read length.
    pub bases2: Option<String>,
    /// The first read's qualities, or the builder's base quality at every base.
    pub quals1: Option<Vec<u8>>,
    /// The second read's qualities, or the builder's base quality at every base.
    pub quals2: Option<Vec<u8>>,
    /// The index of both reads' contig in the header.
    pub contig: usize,
    /// The index of the second read's contig, if not `contig`.
    pub contig2: Option<usize>,
    /// The first read's 1-based start, or 0 for an unmapped read.
    pub start1: usize,
    /// The second read's 1-based start, or 0 for an unmapped read.
    pub start2: usize,
    /// Whether the first read is unmapped.
    pub unmapped1: bool,
    /// Whether the second read is unmapped.
    pub unmapped2: bool,
    /// The first read's CIGAR, or the builder's read length of matches.
    pub cigar1: Option<String>,
    /// The second read's CIGAR, or the builder's read length of matches.
    pub cigar2: Option<String>,
    /// The first read's mapping quality.
    pub mapq1: u8,
    /// The second read's mapping quality.
    pub mapq2: u8,
    /// The first read's strand.
    pub strand1: Strand,
    /// The second read's strand.
    pub strand2: Strand,
    /// Auxiliary fields of both reads.
    pub attrs: Vec<(Tag, Value)>,
}

impl Default for Pair {
    fn default() -> Self {
        Self {
            name: None,
            bases1: None,
            bases2: None,
            quals1: None,
            quals2: None,
            contig: 0,
            contig2: None,
            start1: 0,
            start2: 0,
            unmapped1: false,
            unmapped2: false,
            cigar1: None,
            cigar2: None,
            mapq1: 60,
            mapq2: 60,
            strand1: Strand::Plus,
            strand2: Strand::Minus,
            attrs: Vec::new(),
        }
    }
}

impl Pair {
    /// A default FR pair at two 1-based starts.
    pub fn at(start1: usize, start2: usize) -> Self {
        Self {
            start1,
            start2,
            ..Self::default()
        }
    }

    /// A default FR pair at two 1-based starts with one base repeated along both reads.
    pub fn filled(start1: usize, start2: usize, base: char, read_length: usize) -> Self {
        Self {
            start1,
            start2,
            bases1: Some(base.to_string().repeat(read_length)),
            bases2: Some(base.to_string().repeat(read_length)),
            ..Self::default()
        }
    }
}

/// An unpaired read to build, with the defaults of fgbio's `SamBuilder.addFrag`.
#[derive(Clone, Debug)]
pub struct Frag {
    /// The read's name, or the next sequential name.
    pub name: Option<String>,
    /// The read's bases, or random bases of the builder's read length.
    pub bases: Option<String>,
    /// The read's qualities, or the builder's base quality at every base.
    pub quals: Option<Vec<u8>>,
    /// The index of the read's contig in the header.
    pub contig: usize,
    /// The read's 1-based start, or 0 for an unmapped read.
    pub start: usize,
    /// Whether the read is unmapped.
    pub unmapped: bool,
    /// The read's CIGAR, or the builder's read length of matches.
    pub cigar: Option<String>,
    /// The read's mapping quality.
    pub mapq: u8,
    /// The read's strand.
    pub strand: Strand,
    /// Auxiliary fields of the read.
    pub attrs: Vec<(Tag, Value)>,
}

impl Default for Frag {
    fn default() -> Self {
        Self {
            name: None,
            bases: None,
            quals: None,
            contig: 0,
            start: 0,
            unmapped: false,
            cigar: None,
            mapq: 60,
            strand: Strand::Plus,
            attrs: Vec::new(),
        }
    }
}

impl Frag {
    /// A default forward read at a 1-based start.
    pub fn at(start: usize) -> Self {
        Self {
            start,
            ..Self::default()
        }
    }
}

/// One read of a pair or a fragment, as its builder describes it.
struct Read<'a> {
    name: &'a str,
    bases: &'a str,
    quals: Vec<u8>,
    contig: usize,
    start: usize,
    unmapped: bool,
    cigar: Cigar,
    mapq: u8,
    strand: Strand,
    flags: Flags,
    attrs: &'a [(Tag, Value)],
}

/// Builds alignment records as fgbio's `SamBuilder` does: contigs `chr1` to `chr22`, `chrX`,
/// `chrY`, and `chrM` of 200 Mbp, one read group, sequential names, random bases, and mate
/// information on pairs as htsjdk's `SamPairUtil.setProperPairAndMateInfo` sets it.
///
/// Its header declares coordinate order, the order in which it writes and sources its records.
#[derive(Clone, Debug)]
pub struct SamBuilder {
    read_length: usize,
    base_quality: u8,
    header: sam::Header,
    records: Vec<RecordBuf>,
    counter: usize,
    seed: u64,
}

impl Default for SamBuilder {
    fn default() -> Self {
        Self::new()
    }
}

/// The names of the contigs of fgbio's default sequence dictionary.
pub fn default_contigs() -> Vec<String> {
    (1..=22)
        .map(|index| index.to_string())
        .chain(["X", "Y", "M"].map(String::from))
        .map(|name| format!("chr{name}"))
        .collect()
}

impl SamBuilder {
    /// A builder of reads of length 100 at base quality 30.
    pub fn new() -> Self {
        let length = NonZeroUsize::new(200_000_000).expect("a contig length");
        let mut header = sam::Header::builder().set_header(
            Map::<Header>::builder()
                .set_version(Version::new(1, 6))
                .insert(SORT_ORDER, "coordinate")
                .build()
                .expect("a valid @HD"),
        );
        for name in default_contigs() {
            header = header.add_reference_sequence(name, Map::<ReferenceSequence>::new(length));
        }
        let read_group = Map::<ReadGroup>::builder()
            .insert(SAMPLE, SAMPLE_NAME)
            .build()
            .expect("a valid @RG");
        Self {
            read_length: 100,
            base_quality: 30,
            header: header.add_read_group(READ_GROUP_ID, read_group).build(),
            records: Vec::new(),
            counter: 0,
            seed: 42,
        }
    }

    /// Sets the length of reads given no bases.
    #[must_use]
    pub fn read_length(mut self, read_length: usize) -> Self {
        self.read_length = read_length;
        self
    }

    /// Sets the quality of every base of reads given no qualities.
    #[must_use]
    pub fn base_quality(mut self, base_quality: u8) -> Self {
        self.base_quality = base_quality;
        self
    }

    /// The header of the built records.
    pub fn header(&self) -> &sam::Header {
        &self.header
    }

    /// Every built record, in the order it was added.
    pub fn records(&self) -> &[RecordBuf] {
        &self.records
    }

    /// Adds records built elsewhere, as fgbio's `++=` does.
    pub fn extend(&mut self, records: impl IntoIterator<Item = RecordBuf>) {
        self.records.extend(records);
    }

    /// Adds a pair of reads, the first and then the second, and returns them.
    ///
    /// Each read points at its mate, with the mate's strand, its mapping quality in `MQ`, and its
    /// CIGAR in `MC`; the template length is htsjdk's insert size, and both reads are flagged as
    /// a proper pair when both are mapped to one contig. An unmapped read takes the position of
    /// its mapped mate.
    pub fn add_pair(&mut self, pair: Pair) -> Vec<RecordBuf> {
        let Pair {
            name,
            bases1,
            bases2,
            quals1,
            quals2,
            contig,
            contig2,
            start1,
            start2,
            unmapped1,
            unmapped2,
            cigar1,
            cigar2,
            mapq1,
            mapq2,
            strand1,
            strand2,
            attrs,
        } = pair;
        let name = name.unwrap_or_else(|| self.next_name());
        let bases1 = bases1.unwrap_or_else(|| self.random_bases());
        let bases2 = bases2.unwrap_or_else(|| self.random_bases());
        let mut r1 = record(&Read {
            name: &name,
            quals: self.quals(quals1, &bases1),
            bases: &bases1,
            contig,
            start: start1,
            unmapped: unmapped1 || start1 == 0,
            cigar: self.cigar(cigar1),
            mapq: mapq1,
            strand: strand1,
            flags: Flags::SEGMENTED | Flags::FIRST_SEGMENT,
            attrs: &attrs,
        });
        let mut r2 = record(&Read {
            name: &name,
            quals: self.quals(quals2, &bases2),
            bases: &bases2,
            contig: contig2.unwrap_or(contig),
            start: start2,
            unmapped: unmapped2 || start2 == 0,
            cigar: self.cigar(cigar2),
            mapq: mapq2,
            strand: strand2,
            flags: Flags::SEGMENTED | Flags::LAST_SEGMENT,
            attrs: &attrs,
        });
        set_mate_info(&mut r1, &mut r2);
        let built = vec![r1, r2];
        self.records.extend(built.iter().cloned());
        built
    }

    /// Adds an unpaired read and returns it.
    pub fn add_frag(&mut self, frag: Frag) -> Vec<RecordBuf> {
        let Frag {
            name,
            bases,
            quals,
            contig,
            start,
            unmapped,
            cigar,
            mapq,
            strand,
            attrs,
        } = frag;
        let name = name.unwrap_or_else(|| self.next_name());
        let bases = bases.unwrap_or_else(|| self.random_bases());
        let built = record(&Read {
            name: &name,
            quals: self.quals(quals, &bases),
            bases: &bases,
            contig,
            start,
            unmapped: unmapped || start == 0,
            cigar: self.cigar(cigar),
            mapq,
            strand,
            flags: Flags::empty(),
            attrs: &attrs,
        });
        self.records.push(built.clone());
        vec![built]
    }

    /// The built records in coordinate order, as BAM records.
    ///
    /// Records with no contig come last, and records at one position keep the order they were
    /// added in.
    pub fn to_bam_records(&self) -> Vec<bam::Record> {
        let mut writer = bam::io::Writer::from(Vec::new());
        for record in self.sorted() {
            writer
                .write_alignment_record(&self.header, record)
                .expect("a built record is encoded");
        }
        let bytes = writer.into_inner();
        let mut reader = bam::io::Reader::from(&bytes[..]);
        let mut records = Vec::with_capacity(self.records.len());
        let mut record = bam::Record::default();
        while reader
            .read_record(&mut record)
            .expect("a built record is decoded")
            > 0
        {
            records.push(record.clone());
        }
        records
    }

    /// A source of the built records in coordinate order, for a [`StreamingPileupBuilder`].
    pub fn to_source(&self) -> BuiltRecords {
        let records: Vec<io::Result<bam::Record>> =
            self.to_bam_records().into_iter().map(Ok).collect();
        Records::new(records.into_iter())
    }

    /// A [`StreamingPileupBuilder`] over the built records, piling them up in memory.
    pub fn to_pileup_builder<'f>(&self) -> StreamingPileupBuilder<'f, BuiltRecords> {
        StreamingPileupBuilder::new(self.to_source(), &self.header)
            .expect("the header declares coordinate order")
    }

    /// Writes the records to a BAM in coordinate order, as fgbio's `write` does.
    pub fn write_bam(&self, path: impl AsRef<Path>) -> io::Result<()> {
        let mut writer = bam::io::Writer::new(File::create(path)?);
        writer.write_header(&self.header)?;
        for record in self.sorted() {
            writer.write_alignment_record(&self.header, record)?;
        }
        writer.try_finish()
    }

    /// A copy of a record claiming its mate maps to another contig.
    pub fn with_mate_reference_sequence_id(mut record: RecordBuf, id: usize) -> RecordBuf {
        *record.mate_reference_sequence_id_mut() = Some(id);
        record
    }

    /// A copy of a record without its `MC` tag.
    pub fn without_mate_cigar(mut record: RecordBuf) -> RecordBuf {
        record.data_mut().remove(&Tag::MATE_CIGAR);
        record
    }

    /// A copy of a record with more flag bits set.
    pub fn with_flags(mut record: RecordBuf, bits: u16) -> RecordBuf {
        let flags = record.flags() | Flags::from_bits_truncate(bits);
        *record.flags_mut() = flags;
        record
    }

    /// A copy of a record with other bases, as fgbio's `rec.bases = ...` makes it.
    pub fn with_bases(mut record: RecordBuf, bases: &str) -> RecordBuf {
        *record.sequence_mut() = Sequence::from(bases.as_bytes().to_vec());
        record
    }

    fn sorted(&self) -> Vec<&RecordBuf> {
        let mut records: Vec<&RecordBuf> = self.records.iter().collect();
        records.sort_by_key(|record| {
            (
                record.reference_sequence_id().unwrap_or(usize::MAX),
                record.alignment_start().map_or(usize::MAX, usize::from),
            )
        });
        records
    }

    fn next_name(&mut self) -> String {
        let name = format!("{:04}", self.counter);
        self.counter += 1;
        name
    }

    fn random_bases(&mut self) -> String {
        (0..self.read_length)
            .map(|_| {
                self.seed = self
                    .seed
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                char::from(b"ACGT"[(self.seed >> 62) as usize])
            })
            .collect()
    }

    fn cigar(&self, cigar: Option<String>) -> Cigar {
        let text = cigar.unwrap_or_else(|| format!("{}M", self.read_length));
        sam::record::Cigar::new(text.as_bytes())
            .iter()
            .collect::<Result<Cigar, _>>()
            .unwrap_or_else(|error| panic!("{text} is not a valid CIGAR: {error}"))
    }

    fn quals(&self, quals: Option<Vec<u8>>, bases: &str) -> Vec<u8> {
        quals.unwrap_or_else(|| vec![self.base_quality; bases.len()])
    }
}

/// A record of one read of a pair or a fragment, as fgbio builds it.
fn record(read: &Read<'_>) -> RecordBuf {
    assert_eq!(
        read.bases.len(),
        read.quals.len(),
        "the bases and qualities of {} have different lengths",
        read.name
    );
    assert!(
        read.unmapped || read.bases.len() == read.cigar.read_length(),
        "the bases of {} do not agree with its CIGAR on length",
        read.name
    );
    let mut flags = read.flags;
    if read.unmapped {
        flags |= Flags::UNMAPPED;
    }
    if read.strand == Strand::Minus {
        flags |= Flags::REVERSE_COMPLEMENTED;
    }
    let mut data: Data = read.attrs.iter().cloned().collect();
    data.insert(Tag::READ_GROUP, Value::from(READ_GROUP_ID));
    let mut builder = RecordBuf::builder()
        .set_name(read.name)
        .set_flags(flags)
        .set_sequence(Sequence::from(read.bases.as_bytes().to_vec()))
        .set_quality_scores(QualityScores::from(read.quals.clone()))
        .set_data(data);
    if let Some(start) = Position::new(read.start) {
        builder = builder
            .set_reference_sequence_id(read.contig)
            .set_alignment_start(start);
    }
    let mapq = if read.unmapped { 0 } else { read.mapq };
    if let Some(mapq) = MappingQuality::new(mapq) {
        builder = builder.set_mapping_quality(mapq);
    }
    if !read.unmapped {
        builder = builder.set_cigar(read.cigar.clone());
    }
    builder.build()
}

/// The 5′ end of a mapped record on the reference, for the insert size.
fn five_prime(record: &RecordBuf) -> i64 {
    let position = if record.flags().is_reverse_complemented() {
        record.alignment_end()
    } else {
        record.alignment_start()
    };
    position.map_or(0, |position| usize::from(position) as i64)
}

/// htsjdk's `SamPairUtil.computeInsertSize`.
fn insert_size(first: &RecordBuf, second: &RecordBuf) -> i32 {
    if first.flags().is_unmapped()
        || second.flags().is_unmapped()
        || first.reference_sequence_id() != second.reference_sequence_id()
    {
        return 0;
    }
    let (first, second) = (five_prime(first), five_prime(second));
    let adjustment = if second >= first { 1 } else { -1 };
    (second - first + adjustment) as i32
}

/// The text of a record's CIGAR.
fn cigar_text(record: &RecordBuf) -> String {
    use noodles::sam::alignment::record::cigar::op::Kind;
    let mut text = String::new();
    for op in record.cigar().as_ref() {
        let symbol = match op.kind() {
            Kind::Match => 'M',
            Kind::Insertion => 'I',
            Kind::Deletion => 'D',
            Kind::Skip => 'N',
            Kind::SoftClip => 'S',
            Kind::HardClip => 'H',
            Kind::Pad => 'P',
            Kind::SequenceMatch => '=',
            Kind::SequenceMismatch => 'X',
        };
        text.push_str(&op.len().to_string());
        text.push(symbol);
    }
    text
}

/// Points each read of a pair at the other, as htsjdk's `SamPairUtil.setProperPairAndMateInfo`
/// does with mate CIGARs on and every pair orientation allowed.
fn set_mate_info(r1: &mut RecordBuf, r2: &mut RecordBuf) {
    let (unmapped1, unmapped2) = (r1.flags().is_unmapped(), r2.flags().is_unmapped());
    if unmapped1 && unmapped2 {
        let (reverse1, reverse2) = (
            r1.flags().is_reverse_complemented(),
            r2.flags().is_reverse_complemented(),
        );
        for (record, mate_reverse) in [(&mut *r1, reverse2), (&mut *r2, reverse1)] {
            *record.reference_sequence_id_mut() = None;
            *record.alignment_start_mut() = None;
            *record.mate_reference_sequence_id_mut() = None;
            *record.mate_alignment_start_mut() = None;
            let mut flags = record.flags() | Flags::MATE_UNMAPPED;
            flags.set(Flags::MATE_REVERSE_COMPLEMENTED, mate_reverse);
            flags.remove(Flags::PROPERLY_SEGMENTED);
            *record.flags_mut() = flags;
            record.data_mut().remove(&Tag::MATE_MAPPING_QUALITY);
            record.data_mut().remove(&Tag::MATE_CIGAR);
            *record.template_length_mut() = 0;
        }
    } else if unmapped1 || unmapped2 {
        let (mapped, unmapped) = if unmapped1 { (r2, r1) } else { (r1, r2) };
        *unmapped.reference_sequence_id_mut() = mapped.reference_sequence_id();
        *unmapped.alignment_start_mut() = mapped.alignment_start();

        *mapped.mate_reference_sequence_id_mut() = unmapped.reference_sequence_id();
        *mapped.mate_alignment_start_mut() = unmapped.alignment_start();
        let mut flags = mapped.flags() | Flags::MATE_UNMAPPED;
        flags.set(
            Flags::MATE_REVERSE_COMPLEMENTED,
            unmapped.flags().is_reverse_complemented(),
        );
        flags.remove(Flags::PROPERLY_SEGMENTED);
        *mapped.flags_mut() = flags;
        mapped.data_mut().remove(&Tag::MATE_MAPPING_QUALITY);
        mapped.data_mut().remove(&Tag::MATE_CIGAR);
        *mapped.template_length_mut() = 0;

        *unmapped.mate_reference_sequence_id_mut() = mapped.reference_sequence_id();
        *unmapped.mate_alignment_start_mut() = mapped.alignment_start();
        let mut flags = unmapped.flags();
        flags.set(
            Flags::MATE_REVERSE_COMPLEMENTED,
            mapped.flags().is_reverse_complemented(),
        );
        flags.remove(Flags::MATE_UNMAPPED | Flags::PROPERLY_SEGMENTED);
        *unmapped.flags_mut() = flags;
        let mapq = mapped.mapping_quality().map_or(255, u8::from);
        unmapped
            .data_mut()
            .insert(Tag::MATE_MAPPING_QUALITY, Value::from(mapq));
        unmapped
            .data_mut()
            .insert(Tag::MATE_CIGAR, Value::from(cigar_text(mapped)));
        *unmapped.template_length_mut() = 0;
    } else {
        let proper = r1.reference_sequence_id() == r2.reference_sequence_id();
        let mates = [
            (r2.reference_sequence_id(), r2.alignment_start()),
            (r1.reference_sequence_id(), r1.alignment_start()),
        ];
        let reverse = [
            r2.flags().is_reverse_complemented(),
            r1.flags().is_reverse_complemented(),
        ];
        let mapqs = [
            r2.mapping_quality().map_or(255, u8::from),
            r1.mapping_quality().map_or(255, u8::from),
        ];
        let cigars = [cigar_text(r2), cigar_text(r1)];
        let tlen = insert_size(r1, r2);
        for (index, (record, tlen)) in [(&mut *r1, tlen), (&mut *r2, -tlen)]
            .into_iter()
            .enumerate()
        {
            *record.mate_reference_sequence_id_mut() = mates[index].0;
            *record.mate_alignment_start_mut() = mates[index].1;
            let mut flags = record.flags();
            flags.set(Flags::MATE_REVERSE_COMPLEMENTED, reverse[index]);
            flags.remove(Flags::MATE_UNMAPPED);
            flags.set(Flags::PROPERLY_SEGMENTED, proper);
            *record.flags_mut() = flags;
            record
                .data_mut()
                .insert(Tag::MATE_MAPPING_QUALITY, Value::from(mapqs[index]));
            record
                .data_mut()
                .insert(Tag::MATE_CIGAR, Value::from(cigars[index].clone()));
            *record.template_length_mut() = tlen;
        }
    }
}
