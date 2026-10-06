# streampile

Forward-only pileups streamed from coordinate-sorted BAM records.

A `StreamingPileupBuilder` reads records once, from an indexed file or an unindexed pipe, and piles them up at the 0-based positions asked for, moving forward only.
Each record's CIGAR is decoded once into a footprint, and a pileup borrows the records it holds, so building one allocates nothing per entry, and per-base auxiliary arrays are read in place.
Entries follow htslib: a deletion carries the quality of the read's next base, a read with no stored qualities has quality 255 at every base, and an `N` base is a no-call.
Insertions are reported at both ends of an alignment.

```rust,no_run
use noodles::bam;
use streampile::StreamingPileupBuilder;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut reader = bam::io::reader::Builder::default().build_from_path("reads.bam")?;
    let header = reader.read_header()?;
    let mut builder = StreamingPileupBuilder::new(reader, &header)?
        .min_mapping_quality(20)
        .min_base_quality(30);

    let mut columns = builder.columns("chr1", 0, 1_000)?;
    while let Some(pileup) = columns.next_pileup() {
        let pileup = pileup?;
        for entry in pileup.iter().filter(|entry| entry.passes(30)) {
            let ends = (entry.five_prime_distance(), entry.template_end_distance()?);
            println!("{} {:?} {:?}", pileup.position(), entry.base(), ends);
        }
    }
    builder.close()?;
    Ok(())
}
```

The `template_end_distance()` method counts the template's bases to the 5′ end of the mate of a read in an FR pair, walking the mate's CIGAR from the `MC` tag and never reading the template length (TLEN), and is an error for a read of an FR pair without a usable `MC` tag.
A pair is FR, as `is_fr_pair()` says, when the forward read's aligned 5′ position is at or before the reverse read's, as in htsjdk 5.0.0, so a read of any other pair has no template end.
A pileup's `templates()` sees each template once, calling the bases of overlapping mates into one with the agreement and disagreement strategies of fgbio's `CallOverlappingConsensusBases`, as fgumi implements them.
A `tap` receives every record once, in input order, so records can be written on as they are passed.
The functions `five_prime_distance`, `template_end_distance`, and `is_fr_pair` give the same answers for any noodles alignment record.
A source of records can be any coordinate-sorted stream, and a record type can carry more beside its BAM record.

## Features

- `libdeflate` (default): inflates BGZF blocks with libdeflate, which compiles C code. Build with `default-features = false` to inflate them in pure Rust.
- `testing`: adds `streampile::testing`, whose `SamBuilder` builds test records as fgbio's does and piles them up in memory with a `StreamingPileupBuilder`.

The same core runs the [streampile](https://pypi.org/project/streampile/) Python package, which reads records with pysam.
