# Benchmarks

The `benchmark.py` script piles up every base of a territory and counts, at each base, the reads holding each base at base quality 30 or more, the reads with a deletion there, and the reads with an insertion after it, from reads with mapping quality 20 or more.
It runs each engine in its own process, prints its time and peak resident memory, and prints a digest of its counts.
The `htslib` and `streampile` digests agree unless a read opens with an insertion, which only streampile reports, or has no stored qualities (QUAL `*`), which streampile counts at quality 255 and the `htslib` engine skips; nothing checks them.

```console
uv run python benchmarks/benchmark.py \
    --bam reads.bam --ref reference.fa --intervals territory.bed --threads 4 --runs 3
```

The `--engines` option picks the engines, `--threads` sets the threads pysam decompresses the BAM with, and `--runs` reports each engine's fastest of that many runs.
The same sweep runs on the Rust crate alone, with no Python, printing the same digest, with the threads decompressing the BAM given last:

```console
cargo run --release --example sweep -- reads.bam territory.bed 4
```

The engines are:

- `htslib`: pysam's `AlignmentFile.pileup` over each span, counting each entry in Python.
- `streampile`: `StreamingPileupBuilder.columns` over each span, counting each entry in Python.
- `tabulate`: `streampile.tabulate`, which also groups MNVs and normalizes indels.
- `rebuild`: building each column from scratch from every overlapping read's aligned pairs, on `--rebuild-columns` evenly spaced columns (2,000 by default; 500 in the results below).
- `examples/sweep`: the Rust crate's `StreamingPileupBuilder` over each span, reading the BAM with noodles and counting each entry in Rust.

The BAMs are synthetic, coordinate-sorted, overlapping 2x150 pairs over 600 targets of 205 bases (123 kb), with the per-base tags of a duplex consensus (`cd`, `ce`, `ad`, `bd`, `ae`, `be`, `ac`, `bc`, `aq`, `bq`), about 2.4 kB of tags a read, and occasional mismatches, indels, low base qualities, low mapping qualities, and QC-fail reads.
The 30x BAM holds 54,000 reads (29 MB), and the 800x BAM 1.44 million reads (748 MB, 3.9 GB decompressed).
Results on an Apple M3 Max with 16 cores, the best of three runs, at load averages of 3 to 5 from other work, in seconds:

| Engine           | Columns | 30x, 1 thread | 30x, 4 threads | 800x, 1 thread | 800x, 4 threads | 800x, 4 threads, peak MB |
|------------------|--------:|--------------:|---------------:|---------------:|----------------:|-------------------------:|
| `htslib`         | 123,000 |           2.9 |            2.7 |           86.7 |            83.5 |                       52 |
| `streampile`     | 123,000 |           1.0 |            1.0 |           25.4 |            23.4 |                       57 |
| `tabulate`       | 123,000 |           0.4 |            0.3 |            5.9 |             2.6 |                       39 |
| `rebuild`        |     500 |           0.3 |            0.2 |            6.0 |             4.0 |                       51 |
| `examples/sweep` | 123,000 |           0.2 |            0.1 |            4.3 |             2.6 |                       15 |

At 800x, the `streampile` engine is about 3.5 times as fast as the `htslib` engine, and most of its time is spent counting each of its entries in Python.
Sweeping the same columns through the Python API without touching an entry takes 4.5 seconds with one decompression thread and 2.5 with four, where pysam alone takes 2.7 and 1.0 to iterate the records.
The `tabulate` engine is 2.3 times as fast with four decompression threads, and the Rust sweep 1.7 times, so decompressing the BAM is a large part of each.
Rebuilding each column from aligned pairs costs about 12 ms a column at 800x, or about 25 minutes for this territory.
