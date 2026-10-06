use std::io;

fn invalid_data(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.into())
}

const MATCH: u32 = 0;
const INSERTION: u32 = 1;
const DELETION: u32 = 2;
const SKIP: u32 = 3;
const SOFT_CLIP: u32 = 4;
const HARD_CLIP: u32 = 5;
const PAD: u32 = 6;
const SEQUENCE_MATCH: u32 = 7;
const SEQUENCE_MISMATCH: u32 = 8;

/// What a read holds at one reference position it spans.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Located {
    /// The query offset of the base aligned there.
    Base(u32),
    /// A deletion, with the query offset of the read's next base when one follows.
    Deletion(Option<u32>),
    /// A reference skip, the CIGAR `N` operator.
    Skip,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum BlockKind {
    Aligned,
    Deleted,
    Skipped,
}

/// One reference-consuming CIGAR operator: `query` is the offset of its first base, or for a
/// deletion or skip, of the read's next base.
#[derive(Clone, Copy, Debug)]
struct Block {
    start: i64,
    end: i64,
    query: u32,
    kind: BlockKind,
}

/// The inserted bases a read holds right after the reference position `anchor`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Insertion {
    pub anchor: i64,
    pub offset: u32,
    pub length: u32,
}

/// Where one read sits on the reference, worked out once from its CIGAR.
///
/// The reference-consuming operators tile `start..end` as blocks, and each lookup moves a cursor
/// forward, so a read is located at each position of a forward sweep in constant time, and a
/// long reference skip costs one block, not one entry per skipped position. The buffers keep
/// their capacity when the footprint is filled again for another read.
#[derive(Debug, Default)]
pub(crate) struct Footprint {
    pub start: i64,
    pub end: i64,
    pub query_length: u32,
    blocks: Vec<Block>,
    insertions: Vec<Insertion>,
    block: usize,
    insertion: usize,
}

impl Footprint {
    /// Walks a raw BAM CIGAR from a 0-based start, and returns whether it has a
    /// reference-consuming operator: M, D, N, =, or X.
    ///
    /// A CIGAR whose query length differs from the number of stored bases is refused, as htslib
    /// refuses it, so a deletion is followed by a base exactly when its offset is within the read.
    pub fn fill(&mut self, start: i64, cigar: &[u8], stored_bases: usize) -> io::Result<bool> {
        self.blocks.clear();
        self.insertions.clear();
        self.block = 0;
        self.insertion = 0;
        self.start = start;
        let (ops, rest) = cigar.as_chunks::<4>();
        if !rest.is_empty() {
            return Err(invalid_data("a CIGAR is not a whole number of operators"));
        }
        let mut position = start;
        let mut query: u64 = 0;
        let mut placed = false;
        for op in ops {
            let op = u32::from_le_bytes(*op);
            let length = op >> 4;
            let offset = to_offset(query)?;
            match op & 0xF {
                MATCH | SEQUENCE_MATCH | SEQUENCE_MISMATCH => {
                    placed = true;
                    self.push_block(position, length, offset, BlockKind::Aligned);
                    position += i64::from(length);
                    query += u64::from(length);
                }
                INSERTION => {
                    self.push_insertion(position - 1, offset, length);
                    query += u64::from(length);
                }
                DELETION => {
                    placed = true;
                    self.push_block(position, length, offset, BlockKind::Deleted);
                    position += i64::from(length);
                }
                SKIP => {
                    placed = true;
                    self.push_block(position, length, offset, BlockKind::Skipped);
                    position += i64::from(length);
                }
                SOFT_CLIP => query += u64::from(length),
                HARD_CLIP | PAD => {}
                kind => return Err(invalid_data(format!("invalid CIGAR operation kind {kind}"))),
            }
        }
        self.end = position;
        if stored_bases > 0 && stored_bases as u64 != query {
            return Err(invalid_data("CIGAR and query sequence lengths differ"));
        }
        self.query_length = to_offset(query)?;
        Ok(placed)
    }

    fn push_block(&mut self, start: i64, length: u32, query: u32, kind: BlockKind) {
        if length > 0 {
            self.blocks.push(Block {
                start,
                end: start + i64::from(length),
                query,
                kind,
            });
        }
    }

    fn push_insertion(&mut self, anchor: i64, offset: u32, length: u32) {
        match self.insertions.last_mut() {
            Some(last)
                if last.anchor == anchor
                    && u64::from(last.offset) + u64::from(last.length) == u64::from(offset) =>
            {
                last.length = last.length.saturating_add(length);
            }
            Some(last) if last.anchor == anchor => {
                *last = Insertion {
                    anchor,
                    offset,
                    length,
                }
            }
            _ => self.insertions.push(Insertion {
                anchor,
                offset,
                length,
            }),
        }
    }

    /// What the read holds at a position, moving the cursor forward to it.
    ///
    /// Positions asked of one footprint must never decrease.
    pub fn locate(&mut self, pos: i64) -> Option<Located> {
        while let Some(block) = self.blocks.get(self.block)
            && block.end <= pos
        {
            self.block += 1;
        }
        let block = self.blocks.get(self.block)?;
        if block.start > pos {
            return None;
        }
        Some(match block.kind {
            BlockKind::Aligned => Located::Base(block.query + (pos - block.start) as u32),
            BlockKind::Deleted => {
                Located::Deletion((block.query < self.query_length).then_some(block.query))
            }
            BlockKind::Skipped => Located::Skip,
        })
    }

    /// The insertion right after a position, moving the cursor forward to it.
    ///
    /// Positions asked of one footprint must never decrease.
    pub fn insertion_at(&mut self, pos: i64) -> Option<Insertion> {
        while let Some(insertion) = self.insertions.get(self.insertion)
            && insertion.anchor < pos
        {
            self.insertion += 1;
        }
        self.insertions
            .get(self.insertion)
            .filter(|insertion| insertion.anchor == pos)
            .copied()
    }
}

fn to_offset(query: u64) -> io::Result<u32> {
    u32::try_from(query).map_err(|_| invalid_data("a read is longer than 4,294,967,295 bases"))
}

/// The number of reference bases a raw BAM CIGAR spans.
pub(crate) fn reference_length(cigar: &[u8]) -> i64 {
    let (ops, _) = cigar.as_chunks::<4>();
    ops.iter()
        .map(|op| u32::from_le_bytes(*op))
        .filter(|op| {
            matches!(
                op & 0xF,
                MATCH | DELETION | SKIP | SEQUENCE_MATCH | SEQUENCE_MISMATCH
            )
        })
        .map(|op| i64::from(op >> 4))
        .sum()
}

/// The number of reference bases a SAM CIGAR string spans, such as the value of an `MC` tag.
pub(crate) fn reference_length_of_text(cigar: &[u8]) -> io::Result<i64> {
    let invalid = || invalid_data(format!("invalid CIGAR {}", String::from_utf8_lossy(cigar)));
    let mut length: i64 = 0;
    let mut count: i64 = 0;
    let mut digits = false;
    for &byte in cigar {
        if byte.is_ascii_digit() {
            count = count
                .checked_mul(10)
                .and_then(|n| n.checked_add(i64::from(byte - b'0')))
                .ok_or_else(invalid)?;
            digits = true;
        } else {
            if !digits {
                return Err(invalid());
            }
            match byte {
                b'M' | b'D' | b'N' | b'=' | b'X' => {
                    length = length
                        .checked_add(count)
                        .filter(|&length| length <= i64::from(u32::MAX))
                        .ok_or_else(invalid)?;
                }
                b'I' | b'S' | b'H' | b'P' => {}
                _ => return Err(invalid()),
            }
            count = 0;
            digits = false;
        }
    }
    if digits { Err(invalid()) } else { Ok(length) }
}
