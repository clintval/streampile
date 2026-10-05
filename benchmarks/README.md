# Benchmarks

`benchmark.py` piles up every base of a territory and counts, at each base, the reads holding each base at base quality 30 or more, the reads with a deletion there, and the reads with an insertion after it, from reads with mapping quality 20 or more.
It runs each engine in its own process, prints its time and peak resident memory, and prints a digest of its counts.
The `htslib`, `streampile`, and `rust` digests agree unless a read opens with an insertion, which only streampile reports, or has no stored qualities (QUAL `*`), which streampile counts at quality 255 and the `htslib` engine skips; nothing checks them.

```console
uv run python benchmarks/benchmark.py --bam reads.bam --ref reference.fa --intervals territory.bed --runs 3
```

`--engines` picks the engines, and `--runs` reports each engine's fastest of that many runs.
The same sweep runs without Python, printing the same digest:

```console
cargo run --release --example sweep -- reads.bam territory.bed
```

The engines are:

- `htslib`: pysam's `AlignmentFile.pileup` over each span, counting each entry in Python.
- `streampile`: `StreamingPileupBuilder.columns` over each span, counting each entry in Python.
- `rust`: the Rust `StreamingPileupBuilder` over each span, counting each entry in Rust and handing Python one dictionary per column, through the private `streampile._native` module.
- `tabulate`: `streampile.tabulate`, which also groups MNVs and normalizes indels.
- `rebuild`: building each column from scratch from every overlapping read's aligned pairs, on `--rebuild-columns` evenly spaced columns (2,000 by default; 500 in the results below).

Results on an Apple M3 Max, for a coordinate-sorted, overlap-clipped duplex consensus BAM of 1.43 million reads over a 123 kb territory, about 750 reads deep, the best of three or more runs on a machine busy with other work (load averages of 12 to 70 on 16 cores):

| Engine       | Columns | Seconds | Peak MB |
|--------------|--------:|--------:|--------:|
| `htslib`     | 123,256 |    68.0 |      55 |
| `streampile` | 123,256 |    66.0 |      43 |
| `tabulate`   | 123,256 |    10.0 |      40 |
| `rebuild`    |     500 |     5.8 |      49 |

Both pileup engines spend most of their time in Python, touching each entry: pysam's htslib engine builds columns in C, but counting the same entries in Python costs as much as the builder's per-read lookups.
Counting with `PileupColumn.get_query_sequences`, which stays in C, takes the htslib engine 6.4 seconds.
Rebuilding each column from aligned pairs costs about 12 ms a column, or about 24 minutes for this territory.
`streampile tabulate` writes the whole table, BGZF with a tabix index, in 10.6 seconds at 41 MB peak resident memory.

The real BAM is not on the machine the Rust engine was added on, so the results below are for synthetic BAMs at about the same depth and at 30x.
Each is a coordinate-sorted BAM of overlapping 2x150 pairs over 600 targets of 205 bases (123 kb), with the per-base tags of a duplex consensus (`cd`, `ce`, `ad`, `bd`, `ae`, `be`, `ac`, `bc`, `aq`, `bq`), about 2.4 kB of tags a read, and occasional mismatches, indels, low base qualities, low mapping qualities, and QC-fail reads.
The 30x BAM holds 54,000 reads (29 MB), and the 800x BAM 1.44 million reads (748 MB, 3.9 GB decompressed).
Results on the same Apple M3 Max, the best of three runs, at load averages of 5 to 20 on 16 cores:

| Engine                      | Columns | 30x seconds | 30x peak MB | 800x seconds | 800x peak MB |
|-----------------------------|--------:|------------:|------------:|-------------:|-------------:|
| `htslib`                    | 123,000 |         2.9 |          37 |         87.3 |           51 |
| `streampile`                | 123,000 |         2.8 |          37 |         80.7 |           49 |
| `rust`                      | 123,000 |         0.2 |          37 |          4.0 |           47 |
| `tabulate`                  | 123,000 |         0.7 |          37 |         16.5 |           38 |
| `rebuild`                   |     500 |         0.3 |          49 |          6.0 |           49 |
| `examples/sweep`, no Python | 123,000 |         0.2 |           3 |          4.5 |           14 |

The `rust` engine is about 20 times as fast as the Python `streampile` engine and the `htslib` engine at 800x, with the same digest as both.
The pure Rust sweep is no faster, so calling the engine from Python adds nothing measurable: Python touches each column once, not each entry.
Most of either run is inflating BGZF blocks, which `gzip -dc` alone takes 3.0 seconds to do for the 800x BAM.
