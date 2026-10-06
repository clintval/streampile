# streampile

[![Build Status](https://github.com/clintval/streampile/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/clintval/streampile/actions/workflows/tests.yml?query=branch%3Amain)
[![Python Versions](https://img.shields.io/badge/python-3.11_|_3.12_|_3.13_|_3.14-blue)](https://github.com/clintval/streampile)
[![Language](https://img.shields.io/badge/language-rust-DEA584.svg)](https://www.rust-lang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/clintval/streampile/blob/main/LICENSE)
[![basedpyright](https://img.shields.io/badge/basedpyright-checked-42b983)](https://docs.basedpyright.com/latest/)
[![mypy](https://www.mypy-lang.org/static/mypy_badge.svg)](https://mypy-lang.org/)
[![uv](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/uv/main/assets/badge/v0.json)](https://docs.astral.sh/uv/)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://docs.astral.sh/ruff/)

Forward-only pileups streamed from coordinate-sorted SAM, BAM, and CRAM records, and a table of the alleles at every base.

[pysam](https://github.com/pysam-developers/pysam) reads the records, so streampile reads whatever pysam reads, and Rust piles them up and tabulates them.
The Rust core is also the [`streampile`](https://crates.io/crates/streampile) crate, which reads BAM.

## Installation

```console
pip install streampile
```

## Quickstart

The examples read `tests/data/reads.bam`, whose reads are named for what they carry, such as `plain`, `snv`, `mnv`, `lowqual`, `insertion`, and `deletion`.

### Building Pileups

A `StreamingPileupBuilder` reads records once and piles them up at the 0-based positions you ask for, moving forward only.
Give the `AlignmentFile` threads to decompress the BAM with, which pays at depth, as the [benchmarks](https://github.com/clintval/streampile/blob/main/benchmarks/README.md) show; for CRAM, also give it `reference_filename`.

```pycon
>>> from pysam import AlignmentFile
>>> from streampile import StreamingPileupBuilder
>>>
>>> with (
...     AlignmentFile("tests/data/reads.bam", threads=4) as reads,
...     StreamingPileupBuilder(reads, min_base_quality=30) as builder,
... ):
...     first = builder.pileup("chr1", 10)
...     second = builder.pileup("chr1", 12)
>>>
>>> first.unfiltered_depth, first.filtered_depth, first.bases
(4, 4, ['A', 'T', 'G', 'A'])
>>> second.unfiltered_depth, second.filtered_depth, second.bases
(4, 3, ['G', 'G', 'G'])

```

At 12, the base of `lowqual` is under the quality floor, so it counts toward `unfiltered_depth` only.

Filter reads with `read_filter`, and count each template once with `without_overlaps()`:

```pycon
>>> with (
...     AlignmentFile("tests/data/reads.bam") as reads,
...     StreamingPileupBuilder(reads, read_filter=lambda read: read.is_paired) as builder,
... ):
...     pileup = builder.pileup("chr1", 50)
>>>
>>> pileup.bases, pileup.without_overlaps().bases
(['T', 'T'], ['T'])

```

Each entry also measures its base's distance from the read's 5′ end and to the template's other end, the 5′ end of its mate in an FR pair, which is found from the `MC` tag:

```pycon
>>> with AlignmentFile("tests/data/reads.bam") as reads, StreamingPileupBuilder(reads) as builder:
...     pileup = builder.pileup("chr1", 45)
>>>
>>> [(entry.alignment.flag, entry.five_prime_distance, entry.template_end_distance) for entry in pileup.pileups]
[(0, 5, None), (99, 5, 14), (147, 14, 5)]

```

The unpaired read has no template end, and the two mates of the pair mirror each other.

### Tapping Every Record

`tap` receives every record, in input order, once the builder has moved past it, so a record can be changed, e.g. tagged, while it is piled up.
Each entry's `alignment` is the very record read, so pass `tap=writer.write` to write the changed records to an `AlignmentFile` opened with `template=reads`.
Pileups see each read as it was when the builder read it; changes made after that, including in `read_filter`, reach `tap` and `alignment` but not later pileups.

```pycon
>>> passed = []
>>> with (
...     AlignmentFile("tests/data/reads.bam", threads=4) as reads,
...     StreamingPileupBuilder(reads, tap=passed.append) as builder,
... ):
...     for entry in builder.pileup("chr1", 10).pileups:
...         entry.alignment.set_tag("XB", entry.base)
>>>
>>> len(passed), [(read.query_name, read.get_tag("XB")) for read in passed if read.has_tag("XB")]
(13, [('plain', 'A'), ('snv', 'T'), ('mnv', 'G'), ('lowqual', 'A')])

```

### Sweeping a Territory

`columns` yields the pileup at every position of a span.
Its views, such as `bases`, `qualities`, and the depths, are computed in Rust, so prefer them to looping over `pileups` in Python, which is where the time goes at depth.

```pycon
>>> from collections import Counter
>>>
>>> with (
...     AlignmentFile("tests/data/reads.bam") as reads,
...     StreamingPileupBuilder(reads, min_base_quality=30) as builder,
... ):
...     counts = [(pileup.reference_pos, Counter(pileup.bases)) for pileup in builder.columns("chr1", 9, 13)]
>>>
>>> for pos, count in counts:
...     print(pos, dict(count))
9 {'A': 4}
10 {'A': 2, 'T': 1, 'G': 1}
11 {'C': 3, 'G': 1}
12 {'G': 3}

```

Each pileup is a snapshot that outlives the builder moving on.

### Tabulating Alleles

`tabulate` counts the reads of every allele at every base of a [bedspec](https://github.com/clintval/bedspec) `Territory`.
Alleles are normalized VCF alleles at 1-based positions, so they match a VCF by `CHROM`, `POS`, `REF`, and `ALT`.

```pycon
>>> from bedspec import Bed3
>>> from bedspec import Territory
>>> from pysam import FastaFile
>>> from streampile import tabulate
>>>
>>> territory = Territory([Bed3("chr1", start=9, end=12)])
>>> with (
...     AlignmentFile("tests/data/reads.bam") as reads,
...     FastaFile("tests/data/reference.fa") as reference,
... ):
...     bases = list(tabulate(reads, reference, territory, min_base_quality=30))
>>>
>>> for base in bases:
...     print(base.pos, base.ref, base.depth, base.alts, base.alt_reads, base.alt_fwd, base.alt_rev)
10 A 4 () () () ()
11 A 4 ('T', 'GG') (1, 1) (1, 0) (0, 1)
12 C 4 () () () ()

```

`alt_fwd` and `alt_rev` split each allele's reads by strand, and `no_calls` counts `N` bases apart from the depth.

### Reading a Table

```pycon
>>> from streampile import TabulationReader
>>>
>>> for base in TabulationReader.from_path("tests/data/counts.tsv"):
...     print(base.pos, base.depth, base.alts, base.alt_reads)
21 5 () ()
22 5 () ()
23 4 ('C',) (1,)
24 4 () ()

```

## Command Line

```console
streampile tabulate \
    --bam tests/data/reads.bam \
    --ref tests/data/reference.fa \
    --intervals tests/data/territory.bed \
    --min-base-quality 30 \
    --min-mapping-quality 20 \
    --out counts.tsv
```

```text
##streampile-tabulation=1
##streampile-version=0.1.0
##bam=tests/data/reads.bam
##reference=tests/data/reference.fa
##territory=tests/data/territory.bed
##min_base_quality=30
##min_mapping_quality=20
##exclude_flags=0xf00
#contig	pos	ref	depth	no_calls	ref_reads	ref_fwd	ref_rev	alt_refs	alts	alt_reads	alt_fwd	alt_rev
chr1	21	T	5	0	5	3	2					
chr1	22	G	5	0	5	3	2					
chr1	23	C	4	0	3	2	1	CA	C	1	0	1
chr1	24	A	4	0	3	2	1					
```

`--threads` sets the threads decompressing the BAM and compressing a BGZF table.
Write to a `.gz` path with `--index tbi` to compress and index the table, then read a region back:

```python
for base in TabulationReader.query("counts.tsv.gz", "chr1", 20, 24):
    print(base.pos, base.depth, base.alts)
```

## Development and Testing

See the [contributing guide](https://github.com/clintval/streampile/blob/main/CONTRIBUTING.md) for more information.

The streaming design follows the `StreamingPileupBuilder` of [fgbio](https://github.com/fulcrumgenomics/fgbio) of which I was also the author.
See [NOTICE](https://github.com/clintval/streampile/blob/main/NOTICE) for official attribution.
