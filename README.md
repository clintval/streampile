# streampile

[![Build Status](https://github.com/clintval/streampile/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/clintval/streampile/actions/workflows/tests.yml?query=branch%3Amain)
[![Python Versions](https://img.shields.io/badge/python-3.11_|_3.12_|_3.13_|_3.14-blue)](https://github.com/clintval/streampile)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/clintval/streampile/blob/main/LICENSE)
[![basedpyright](https://img.shields.io/badge/basedpyright-checked-42b983)](https://docs.basedpyright.com/latest/)
[![mypy](https://www.mypy-lang.org/static/mypy_badge.svg)](https://mypy-lang.org/)
[![uv](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/uv/main/assets/badge/v0.json)](https://docs.astral.sh/uv/)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://docs.astral.sh/ruff/)

Forward-only pileups streamed from coordinate-sorted SAM, BAM, and CRAM records.

## Installation

The package can be installed with `pip`:

```console
pip install streampile
```

## Quickstart

### Building Pileups

A `StreamingPileupBuilder` reads coordinate-sorted records once, from start to finish, and piles them up at the positions you ask for.
A position may repeat or move forward, but never back.
Positions are 0-based, as in pysam.

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
>>> first.filtered_depth, first.get_query_sequences
(4, ['A', 'T', 'G', 'A'])
>>> second.filtered_depth, second.get_query_sequences
(3, ['G', 'G', 'G'])

```

Each entry of `pileup.pileups` holds a read and its base, deletion, reference skip (`N`), or insertion at the position.
A skip holds no base or quality: it counts in `unfiltered_depth`, as in htslib, but never in `filtered_depth` or the bases and qualities of a pileup.
Each read's CIGAR is walked once, when the builder first reaches it, so a pileup costs one lookup per read.
Pass `tap`, e.g. `tap=writer.write`, to be handed every record, in input order, once the builder has moved past it.
Pass `read_filter`, e.g. `read_filter=lambda read: read.is_proper_pair`, to leave more reads out of pileups, after the built-in filters; a read it rejects still goes to `tap`.

### Sweeping a Territory

`columns` yields the pileup at every position of a span, covered or not.

```pycon
>>> with AlignmentFile("tests/data/reads.bam") as reads, StreamingPileupBuilder(reads) as builder:
...     sum(pileup.unfiltered_depth for pileup in builder.columns("chr1", 0, 20))
73

```

### Tabulating Alleles

`tabulate` counts the reads of every allele at every base of a territory, one `TabulatedBase` per base, covered or not.
Within a read, adjacent mismatches are grouped into one MNV, and indels are anchored, trimmed, and left-aligned as normalized VCF alleles.
Positions are 1-based, as in VCF, so alleles can be matched to a VCF by `CHROM`, `POS`, `REF`, and `ALT`.

```pycon
>>> from pysam import FastaFile
>>> from streampile import tabulate
>>>
>>> with (
...     AlignmentFile("tests/data/reads.bam") as reads,
...     FastaFile("tests/data/reference.fa") as reference,
... ):
...     bases = list(tabulate(reads, reference, [("chr1", 9, 12)], min_base_quality=30))
>>>
>>> for base in bases:
...     print(base.pos, base.ref, base.depth, base.ref_reads, base.alt_refs, base.alts, base.alt_reads)
10 A 4 4 () () ()
11 A 4 2 ('A', 'AC') ('T', 'GG') (1, 1)
12 C 4 3 () () ()

```

The alleles anchored at a base sit in parallel tuples, as VCF pairs `ALT` with `AD`: allele `i` is `alt_refs[i]` to `alts[i]`, seen in `alt_reads[i]` reads.
A read counts toward `depth` at a base when it holds an aligned base there at the quality floor, or when the base lies inside an allele it was counted for, such as the deleted bases of a deletion or the second base of an MNV.
It is counted once: as a reference read, for the allele it has anchored at the base, or in `depth` alone.
Overlapping mates are both counted, so clip overlaps first to count each molecule once.

### Reading a Table

`TabulationReader` reads a table written by `streampile tabulate` or `TabulationWriter`, compressed or not, back into typed records.

```pycon
>>> from streampile import TabulationReader
>>>
>>> for base in TabulationReader.from_path("tests/data/counts.tsv"):
...     print(base.pos, base.depth, base.alt_refs, base.alts, base.alt_reads)
21 5 () () ()
22 5 () () ()
23 4 ('CA',) ('C',) (1,)
24 4 () () ()

```

## Command Line

`streampile tabulate` writes the same bases to a table:

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
contig	pos	ref	depth	ref_reads	alt_refs	alts	alt_reads
chr1	21	T	5	5			
chr1	22	G	5	5			
chr1	23	C	4	3	CA	C	1
chr1	24	A	4	3			
```

A table whose path ends in `.gz` or `.bgz` is written as BGZF with [pybgzf](https://github.com/clintval/pybgzf), and `--index tbi` or `--index csi` indexes it by contig and position as it is written.

## Development and Testing

See the [contributing guide](https://github.com/clintval/streampile/blob/main/CONTRIBUTING.md) for more information.

The streaming design follows the `StreamingPileupBuilder` of [fgbio](https://github.com/fulcrumgenomics/fgbio); see [NOTICE](https://github.com/clintval/streampile/blob/main/NOTICE).
