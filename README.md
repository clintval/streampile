# streampile

[![Build Status](https://github.com/clintval/streampile/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/clintval/streampile/actions/workflows/tests.yml?query=branch%3Amain)
[![Python Versions](https://img.shields.io/badge/python-3.11_|_3.12_|_3.13_|_3.14-blue)](https://github.com/clintval/streampile)
[![Language](https://img.shields.io/badge/language-rust-DEA584.svg)](https://www.rust-lang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/clintval/streampile/blob/main/LICENSE)
[![basedpyright](https://img.shields.io/badge/basedpyright-checked-42b983)](https://docs.basedpyright.com/latest/)
[![mypy](https://www.mypy-lang.org/static/mypy_badge.svg)](https://mypy-lang.org/)
[![uv](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/uv/main/assets/badge/v0.json)](https://docs.astral.sh/uv/)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://docs.astral.sh/ruff/)

Forward-only pileups streamed from coordinate-sorted BAM and CRAM records, and a table of the alleles at every base.

## Installation

```console
pip install streampile
```

The streaming pileup core is also a Rust crate, which the Python package will build on:

```console
cargo add streampile
```

## Quickstart

### Building Pileups

A `StreamingPileupBuilder` reads records once and piles them up at the 0-based positions you ask for, moving forward only.

```pycon
>>> from pysam import AlignmentFile
>>> from streampile import StreamingPileupBuilder
>>>
>>> with (
...     AlignmentFile("tests/data/reads.bam") as reads,
...     StreamingPileupBuilder(reads, min_base_quality=30) as builder,
... ):
...     first = builder.pileup("chr1", 10)
...     second = builder.pileup("chr1", 12)
>>>
>>> first.filtered_depth, first.bases
(4, ['A', 'T', 'G', 'A'])
>>> second.filtered_depth, second.bases
(3, ['G', 'G', 'G'])

```

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

Pass `tap=writer.write` to receive every record, in input order, once the builder has moved past it.

### Sweeping a Territory

`columns` yields the pileup at every position of a span.

```pycon
>>> with AlignmentFile("tests/data/reads.bam") as reads, StreamingPileupBuilder(reads) as builder:
...     sum(pileup.unfiltered_depth for pileup in builder.columns("chr1", 0, 20))
73

```

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
...     print(base.pos, base.ref, base.depth, base.alts, base.alt_reads)
10 A 4 () ()
11 A 4 ('T', 'GG') (1, 1)
12 C 4 () ()

```

Each base also splits its reads by strand and counts no-calls apart from its depth.

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

Write to a `.gz` path with `--index tbi` to compress and index the table, then read a region back:

```python
for base in TabulationReader.query("counts.tsv.gz", "chr1", 20, 24):
    print(base.pos, base.depth, base.alts)
```

## Development and Testing

See the [contributing guide](https://github.com/clintval/streampile/blob/main/CONTRIBUTING.md) for more information.

The streaming design follows the `StreamingPileupBuilder` of [fgbio](https://github.com/fulcrumgenomics/fgbio) of which I was also the author.
See [NOTICE](https://github.com/clintval/streampile/blob/main/NOTICE) for official attribution.
