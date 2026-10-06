//! Sweeps every base of a BED territory of a BAM, as `benchmarks/benchmark.py` does, and prints
//! the time, the number of columns, and the digest of the counts the Python engines print, with
//! as many threads decompressing the BAM as given, or one.
//!
//! ```console
//! cargo run --release --example sweep -- reads.bam territory.bed 4
//! ```

use std::collections::BTreeMap;
use std::error::Error;
use std::fmt::Write as _;
use std::fs::{self, File};
use std::io::Read;
use std::num::NonZero;
use std::time::Instant;

use noodles::sam::alignment::record::Flags;
use noodles::{bam, bgzf};
use sha2::{Digest, Sha256};
use streampile::{EntryKind, StreamingPileupBuilder};

const MIN_MAPPING_QUALITY: u8 = 20;
const MIN_BASE_QUALITY: u8 = 30;
const EXCLUDE_FLAGS: u16 = 0xF04;

fn main() -> Result<(), Box<dyn Error>> {
    let mut args = std::env::args().skip(1);
    let (Some(bam_path), Some(bed_path)) = (args.next(), args.next()) else {
        return Err("usage: sweep <reads.bam> <territory.bed> [threads]".into());
    };
    let threads = args
        .next()
        .map_or(Ok(1), |threads| threads.parse::<usize>())?;
    let started = Instant::now();
    let file = File::open(&bam_path)?;
    let decoder: Box<dyn Read> = match NonZero::new(threads).filter(|threads| threads.get() > 1) {
        Some(workers) => Box::new(bgzf::io::MultithreadedReader::with_worker_count(
            workers, file,
        )),
        None => Box::new(bgzf::io::Reader::new(file)),
    };
    let mut reader = bam::io::Reader::from(decoder);
    let header = reader.read_header()?;
    let spans = territory(&fs::read_to_string(&bed_path)?, &header)?;

    let mut builder = StreamingPileupBuilder::new(reader, &header)?
        .min_mapping_quality(MIN_MAPPING_QUALITY)
        .exclude_flags(Flags::from_bits_retain(EXCLUDE_FLAGS));
    let mut digest = Sha256::new();
    let mut columns = 0_usize;
    for (reference_sequence_id, start, end) in spans {
        let mut sweep = builder.columns_at(reference_sequence_id, start, end)?;
        while let Some(pileup) = sweep.next_pileup() {
            let mut counts: BTreeMap<u8, usize> = BTreeMap::new();
            for entry in pileup?.iter() {
                if entry.kind() == EntryKind::Insertion {
                    *counts.entry(b'+').or_default() += 1;
                } else if entry.passes(MIN_BASE_QUALITY) {
                    let key = if entry.is_deletion() {
                        b'-'
                    } else {
                        entry.base().unwrap_or(b'N')
                    };
                    *counts.entry(key).or_default() += 1;
                }
            }
            let items: Vec<String> = counts
                .iter()
                .map(|(key, count)| format!("('{}', {count})", char::from(*key)))
                .collect();
            digest.update(format!("[{}]", items.join(", ")).as_bytes());
            columns += 1;
        }
    }
    builder.close()?;
    let seconds = started.elapsed().as_secs_f64();
    let digest = digest
        .finalize()
        .iter()
        .fold(String::new(), |mut hex, byte| {
            let _ = write!(hex, "{byte:02x}");
            hex
        });
    println!("rust-example\t{seconds:.1}\t{columns}\t{}", &digest[..12]);
    Ok(())
}

/// A contig's index in the header, and a 0-based start and end on it.
type Span = (usize, usize, usize);

/// The spans of a BED file in the order of the header's contigs, overlapping or abutting spans
/// joined, as a `bedspec.Territory` joins them.
fn territory(bed: &str, header: &noodles::sam::Header) -> Result<Vec<Span>, Box<dyn Error>> {
    let mut spans = Vec::new();
    for line in bed
        .lines()
        .filter(|line| !(line.is_empty() || line.starts_with(['#', 't', 'b'])))
    {
        let mut fields = line.split('\t');
        let (Some(contig), Some(start), Some(end)) = (fields.next(), fields.next(), fields.next())
        else {
            return Err(format!("not a BED line: {line}").into());
        };
        let id = header
            .reference_sequences()
            .get_index_of(contig.as_bytes())
            .ok_or_else(|| format!("contig {contig} is not in the header"))?;
        spans.push((id, start.parse::<usize>()?, end.parse::<usize>()?));
    }
    spans.sort_unstable();
    let mut joined: Vec<Span> = Vec::new();
    for (id, start, end) in spans {
        match joined.last_mut() {
            Some(last) if last.0 == id && start <= last.2 => last.2 = last.2.max(end),
            _ if start < end => joined.push((id, start, end)),
            _ => {}
        }
    }
    Ok(joined)
}
