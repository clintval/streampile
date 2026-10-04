# Benchmarks

`benchmark.py` piles up every base of a territory and counts, at each base, the reads holding each base at base quality 30 or more, the reads with a deletion there, and the reads with an insertion after it, from reads with mapping quality 20 or more.
It runs each engine in its own process, prints its time and peak resident memory, and prints a digest of its counts, which must agree between engines.

```console
uv run python benchmarks/benchmark.py --bam reads.bam --ref reference.fa --intervals territory.bed
```

The engines are:

- `htslib`: pysam's `AlignmentFile.pileup` over each span, counting each entry in Python.
- `streampile`: `StreamingPileupBuilder.columns` over each span, counting each entry in Python.
- `tabulate`: `streampile.tabulate`, which also groups MNVs and normalizes indels.
- `rebuild`: building each column from scratch from every overlapping read's aligned pairs, on `--rebuild-columns` evenly spaced columns (500 by default).

Results on an Apple M3 Max, for a coordinate-sorted, overlap-clipped duplex consensus BAM of 1.43 million reads over a 123 kb territory, about 750 reads deep:

| Engine       | Columns | Seconds | Peak MB |
|--------------|--------:|--------:|--------:|
| `htslib`     | 123,256 |    63.1 |      54 |
| `streampile` | 123,256 |    60.3 |      47 |
| `tabulate`   | 123,256 |    10.4 |      38 |
| `rebuild`    |     500 |     5.6 |      48 |

Both pileup engines spend most of their time in Python, touching each entry: pysam's htslib engine builds columns in C, but counting the same entries in Python costs as much as the builder's per-read lookups.
Counting with `PileupColumn.get_query_sequences`, which stays in C, takes the htslib engine 6.1 seconds.
Rebuilding each column from aligned pairs costs about 11 ms a column, or about 23 minutes for this territory.
