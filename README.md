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
It filters reads as `tabulate` does, by `min_mapping_quality` and `exclude_flags`, and by default leaves out secondary, supplementary, duplicate, and QC-fail reads, where htslib keeps supplementary reads and fgbio keeps QC-fail reads.
Its filtered views leave out bases under quality 13, as pysam's `pileup()` and `samtools mpileup` do.

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

Each entry of `pileup.pileups` holds a read and its base, deletion, reference skip (the CIGAR `N` operator), or insertion at the position; `is_no_call` marks a base entry holding an `N` base.
A skip holds no base or quality: it counts in `unfiltered_depth`, as in htslib, but never in `filtered_depth`, `bases`, or `qualities`.
`bases` and `qualities` list only the bases at the quality floor, in the same order, where pysam's `get_query_sequences()` and `get_query_qualities()` list every entry.
An insertion is reported at the position before it, at either end of an alignment: one that opens it, before the first aligned base, and one that closes it, at the last; htslib reports only the closing one, and fgbio only the opening one.
Each read's CIGAR is walked once, when the builder first reaches it, so a pileup costs one lookup per read.
Pass `tap`, e.g. `tap=writer.write`, to be handed every record, in input order, once the builder has moved past it.
Keeping input order holds every read behind the longest read still in a pileup, so a `tap` costs memory with long or spliced reads; without one, a read is dropped once passed.
Pass `read_filter`, e.g. `read_filter=lambda read: read.is_proper_pair`, to leave more reads out of pileups, after the built-in filters; a read it rejects still goes to `tap`.
Mapped pairs, which fgbio keeps by default, are kept with `read_filter=lambda read: read.is_paired and not read.mate_is_unmapped`; no filter of the positions outside an FR pair's insert is provided.
Overlapping mates are both piled up; `pileup.without_overlaps()` keeps one read per template: the first in the input whose base or deletion there is at the quality floor, or else the first, so a mate's skip or low-quality base never hides the other mate's base.

### Sweeping a Territory

`columns` yields the pileup at every position of a span, covered or not.

```pycon
>>> with AlignmentFile("tests/data/reads.bam") as reads, StreamingPileupBuilder(reads) as builder:
...     sum(pileup.unfiltered_depth for pileup in builder.columns("chr1", 0, 20))
73

```

### Tabulating Alleles

`tabulate` counts the reads of every allele at every base of a territory, one `TabulatedBase` per base, covered or not.
The territory is a [bedspec](https://github.com/clintval/bedspec) `Territory`, e.g. `Territory(BedReader.from_path[Bed3N]("territory.bed"))`, whose overlapping and abutting spans are joined, and its bases come in the order of the alignment header's contigs.
Within a read, adjacent mismatches are grouped into one MNV, and indels are anchored, trimmed, and left-aligned as normalized VCF alleles.
Positions are 1-based, as in VCF, so alleles can be matched to a VCF by `CHROM`, `POS`, `REF`, and `ALT`.
Like the builder, it leaves out secondary, supplementary, duplicate, and QC-fail reads by default, but it counts bases of any quality unless given a `min_base_quality`, where the builder's floor is 13.

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
...     print(base.pos, base.ref, base.depth, base.ref_reads, base.alt_refs, base.alts, base.alt_reads)
10 A 4 4 () () ()
11 A 4 2 ('A', 'AC') ('T', 'GG') (1, 1)
12 C 4 3 () () ()

```

The alleles anchored at a base sit in parallel tuples, as VCF pairs `ALT` with `AD`: allele `i` is `alt_refs[i]` to `alts[i]`, seen in `alt_reads[i]` reads.
Reads are also split by the strand they are mapped to: `ref_reads` is `ref_fwd + ref_rev`, and `alt_reads[i]` is `alt_fwd[i] + alt_rev[i]`.
A read counts toward `depth` at a base when it holds an aligned base there at the quality floor, other than an `N`, over a reference base of A, C, G, or T, or when the base lies inside an allele it was counted for, such as the deleted bases of a deletion or the second base of an MNV.
It is counted once: as a reference read, for the allele it has anchored at the base, or in `depth` alone.
As in the builder's `filtered_depth`, and in pysam's `pileup()` and `samtools mpileup`, a deletion is judged by the quality of the read's next base.
A read that skips over a base with the CIGAR `N` operator observed no base there, so it is not counted at all, unlike in a pileup's `unfiltered_depth`.
A read holding an `N` base, a no-call, is not informative either: it is counted in `no_calls`, not `depth`, so the molecular depth at a base is `depth + no_calls`, less any reads left out by the quality floor or for an allele they could not be counted for.
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

A table whose path ends in `.gz`, `.bgz`, or `.bgzf` is written as BGZF with [pybgzf](https://github.com/clintval/pybgzf), and `--index tbi` or `--index csi` indexes it by contig and position as it is written.
`TabulationReader.query` reads the rows of a region of an indexed table, 0-based and half-open as in BED, back into typed records:

```python
from streampile import TabulationReader

for base in TabulationReader.query("counts.tsv.gz", "chr1", 20, 24):
    print(base.pos, base.depth, base.no_calls, base.alts, base.alt_reads)
```

A table describes itself.
It opens with `##key=value` lines: its format version, `streampile-tabulation`, the streampile version that wrote it, and every parameter it was tabulated with.
A header line starting `#contig` follows, as VCF's starts `#CHROM`, so an index skips every line before the rows by their `#`.
`TabulationReader` holds the `##` lines in its `metadata` and refuses a table of another format version.
Within a format version, columns are only ever appended, each keeping its meaning: a reader of an earlier streampile reads a later table, keeping the columns it does not know, as text, in `TabulatedBase.extra`.
Any other change to the columns is a new format version.

## Development and Testing

See the [contributing guide](https://github.com/clintval/streampile/blob/main/CONTRIBUTING.md) for more information.

The streaming design follows the `StreamingPileupBuilder` of [fgbio](https://github.com/fulcrumgenomics/fgbio); see [NOTICE](https://github.com/clintval/streampile/blob/main/NOTICE).
