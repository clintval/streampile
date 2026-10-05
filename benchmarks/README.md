# Benchmarks

`benchmark.py` piles up every base of a territory and counts, at each base, the reads holding each base at base quality 30 or more, the reads with a deletion there, and the reads with an insertion after it, from reads with mapping quality 20 or more.
It runs each engine in its own process, prints its time and peak resident memory, and prints a digest of its counts.
The `htslib` and `streampile` digests agree unless a read opens with an insertion, which only streampile reports, or has no stored qualities (QUAL `*`), which streampile counts at quality 255 and the `htslib` engine skips; nothing checks them.

```console
uv run python benchmarks/benchmark.py --bam reads.bam --ref reference.fa --intervals territory.bed
```

The engines are:

- `htslib`: pysam's `AlignmentFile.pileup` over each span, counting each entry in Python.
- `streampile`: `StreamingPileupBuilder.columns` over each span, counting each entry in Python.
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
