# streampile

[![Build Status](https://github.com/clintval/streampile/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/clintval/streampile/actions/workflows/tests.yml?query=branch%3Amain)
[![Python Versions](https://img.shields.io/badge/python-3.11_|_3.12_|_3.13_|_3.14-blue)](https://github.com/clintval/streampile)
[![Language](https://img.shields.io/badge/language-rust-DEA584.svg)](https://www.rust-lang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/clintval/streampile/blob/main/LICENSE)
[![basedpyright](https://img.shields.io/badge/basedpyright-checked-42b983)](https://docs.basedpyright.com/latest/)
[![mypy](https://www.mypy-lang.org/static/mypy_badge.svg)](https://mypy-lang.org/)
[![uv](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/uv/main/assets/badge/v0.json)](https://docs.astral.sh/uv/)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://docs.astral.sh/ruff/)

Forward-only pileups of coordinate-sorted SAM, BAM, and CRAM, and the alleles at every base.

## Installation

```console
pip install streampile
```

## Quickstart

### Building Pileups

A `StreamingPileupBuilder` reads records once and piles them up at the 0-based positions you ask for, moving forward only.
Give the `AlignmentFile` threads to decompress the BAM with, which pays at depth; for CRAM, also give it `reference_filename`.

```python
from pysam import AlignmentFile
from streampile import StreamingPileupBuilder

with (
    AlignmentFile("tests/data/reads.bam", threads=4) as reads,
    StreamingPileupBuilder(reads, min_base_quality=30) as builder,
):
    first = builder.pileup("chr1", 10)
    second = builder.pileup("chr1", 12)

assert first.unfiltered_depth == 4
assert first.filtered_depth == 4
assert first.bases == ["A", "T", "G", "A"]

assert second.unfiltered_depth == 4
assert second.filtered_depth == 3
assert second.bases == ["G", "G", "G"]
```

At position 12, the base of `lowqual` is under the quality floor, so it counts toward `unfiltered_depth` only.

Filter reads with `read_filter`, and see each template once with `templates()`, which calls the bases of overlapping mates into one as fgbio's `CallOverlappingConsensusBases` does.
Overlapping mates read the same bases of their template twice: they agree where they read the same base, and disagree where they read different ones.
The mates of `pair` overlap from position 45 to 54, agree on `T` at position 50, and disagree at position 52, where the first reads `G` at Q40 and the second `A` at Q20:

```text
position  40   45   50   55
pair  99  GATCGGATCCTAGGC
pair 147       GATCCTAAGCTAGCT
                    ^ ^
```

By default, an agreeing base is called at the sum of the two qualities, and a disagreeing one as the base of the higher quality at their difference:

```python
with (
    AlignmentFile("tests/data/reads.bam") as reads,
    StreamingPileupBuilder(reads, read_filter=lambda read: read.is_paired) as builder,
):
    agreeing = builder.pileup("chr1", 50)
    disagreeing = builder.pileup("chr1", 52)

agreed = agreeing.templates()[0]
disagreed = disagreeing.templates()[0]

assert (agreeing.bases, agreeing.qualities) == (["T", "T"], [40, 40])
assert (agreed.base, agreed.qual) == ("T", 80)

assert (disagreeing.bases, disagreeing.qualities) == (["G", "A"], [40, 20])
assert (disagreed.base, disagreed.qual) == ("G", 20)
```

Other strategies mask a disagreement, or keep the base of the higher quality as it is:

```python
from streampile import DisagreementStrategy

masked = disagreeing.templates(disagreement=DisagreementStrategy.mask_both)[0]
kept = disagreeing.templates(disagreement=DisagreementStrategy.mask_lower_qual)[0]

assert (masked.base, masked.qual) == ("N", 2)
assert (kept.base, kept.qual) == ("G", 40)
```

A mate whose base is under the quality floor does not vote, as it is left out of `bases`, so it neither masks nor lowers the other mate's base:

```python
with (
    AlignmentFile("tests/data/reads.bam") as reads,
    StreamingPileupBuilder(reads, min_base_quality=30) as builder,
):
    floored = builder.pileup("chr1", 52)

template = floored.templates(disagreement=DisagreementStrategy.mask_both)[0]

assert floored.bases == ["G"]
assert (template.base, template.qual) == ("G", 40)
```

A template also holds its `reads`, its strand, and its distances to both of its ends:

```python
template = disagreeing.templates()[0]

assert not template.is_reverse
assert template.five_prime_distance == 12
assert template.template_end_distance == 7
```

Each entry also counts the template's bases from its base to the read's 5′ end and to the template's other end, the 5′ end of its mate in an FR pair, whose CIGAR comes from the `MC` tag:

```python
with (
    AlignmentFile("tests/data/reads.bam") as reads,
    StreamingPileupBuilder(reads) as builder,
):
    pileup = builder.pileup("chr1", 45)

distances = [
    (entry.alignment.flag, entry.five_prime_distance, entry.template_end_distance)
    for entry in pileup.pileups
]

assert distances == [(0, 5, None), (99, 5, 14), (147, 14, 5)]
```

The unpaired read has no template end, and the two mates of the pair mirror each other.
A read has one only when `is_fr_pair` calls its pair FR, as htsjdk 5.0.0 does: the forward read's aligned 5′ position is at or before the reverse read's.

### Tapping Every Record

The function given as `tap` receives every record, in input order, once the builder has moved past it, so a record can be changed, e.g. tagged, while it is piled up.
Each entry's `alignment` is the very record read, so pass `tap=writer.write` to write the changed records to an `AlignmentFile` opened with `template=reads`.
Pileups see each read as it was when the builder read it; changes made after that, including in `read_filter`, reach `tap` and `alignment` but not later pileups.
Here every read is written to a new BAM, and a read whose base at position 10 is not the reference `A` is tagged `XV` with that base.

```python
from tempfile import mkdtemp

tagged = f"{mkdtemp()}/tagged.bam"

with (
    AlignmentFile("tests/data/reads.bam") as reads,
    AlignmentFile(tagged, "wb", template=reads) as writer,
    StreamingPileupBuilder(reads, tap=writer.write) as builder,
):
    for entry in builder.pileup("chr1", 10).pileups:
        if entry.base not in (None, "A"):
            entry.alignment.set_tag("XV", entry.base)

with AlignmentFile(tagged) as written:
    variants = {read.query_name: read.get_tag("XV") for read in written if read.has_tag("XV")}

assert variants == {"snv": "T", "mnv": "G"}
```

### Sweeping a Territory

The `columns` method yields the pileup at every position of a span.
Its views, such as `bases`, `qualities`, and the depths, are computed in Rust, so prefer them to looping over `pileups` in Python, which is where the time goes at depth.

```python
from collections import Counter

with (
    AlignmentFile("tests/data/reads.bam") as reads,
    StreamingPileupBuilder(reads, min_base_quality=30) as builder,
):
    counts = {
        pileup.reference_pos: Counter(pileup.bases) for pileup in builder.columns("chr1", 9, 13)
    }

assert counts == {
    9: {"A": 4},
    10: {"A": 2, "T": 1, "G": 1},
    11: {"C": 3, "G": 1},
    12: {"G": 3},
}
```

Each pileup is a snapshot that outlives the builder moving on.

### Tabulating Alleles

The `tabulate` function counts the reads of every allele at every base of a [bedspec](https://github.com/clintval/bedspec) `Territory`.
Alleles are normalized VCF alleles at 1-based positions, so they match a VCF by `CHROM`, `POS`, `REF`, and `ALT`.

```python
from bedspec import Bed3
from bedspec import Territory
from pysam import FastaFile
from streampile import tabulate

territory = Territory([Bed3("chr1", start=9, end=12)])

with (
    AlignmentFile("tests/data/reads.bam") as reads,
    FastaFile("tests/data/reference.fa") as reference,
):
    bases = list(tabulate(reads, reference, territory, min_base_quality=30))

alleles = [
    (base.pos, base.ref, base.alts, base.alt_reads, base.alt_fwd, base.alt_rev) for base in bases
]

assert alleles == [
    (10, "A", (), (), (), ()),
    (11, "A", ("T", "GG"), (1, 1), (1, 0), (0, 1)),
    (12, "C", (), (), (), ()),
]
```

The `alt_fwd` and `alt_rev` fields split each allele's reads by strand, and `no_calls` counts `N` bases apart from the depth.

### Writing a Table

A path ending in `.gz` is written with BGZF, indexed as it is written when given an index, and `metadata` adds `##key=value` lines above the header:

```python
from pybgzf import IndexFormat
from streampile import TabulationWriter

table = f"{mkdtemp()}/counts.tsv.gz"

with TabulationWriter.from_path(
    table, index=IndexFormat.TBI, metadata={"min_base_quality": 30}
) as writer:
    writer.write_all(bases)
```

### Reading a Table

```python
from streampile import TabulationReader

rows = [
    (base.pos, base.depth, base.alts, base.alt_reads) for base in TabulationReader.from_path(table)
]

assert rows == [(10, 4, (), ()), (11, 4, ("T", "GG"), (1, 1)), (12, 4, (), ())]
```

With the index, a region can be read alone; it is half-open and 0-based, as in BED:

```python
region = TabulationReader.query(table, "chr1", 10, 11)
rows = [(base.pos, base.depth, base.alts, base.alt_reads) for base in region]

assert rows == [(11, 4, ("T", "GG"), (1, 1))]
```

## Command Line

```console
streampile tabulate \
    --bam tests/data/reads.bam \
    --ref tests/data/reference.fa \
    --intervals tests/data/territory.bed \
    --min-base-quality 30 \
    --min-mapping-quality 20 \
    --out counts.tsv.gz \
    --index tbi
```

This writes the table compressed with BGZF and its tabix index, `counts.tsv.gz.tbi`; decompressed, it reads:

```text
##streampile-tabulation=1
##streampile-version=0.2.0
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

The `--threads` option sets the threads decompressing the BAM and compressing a BGZF table: by default, the CPUs available, at most 8.

The examples read `tests/data/reads.bam`, whose reads are named for what they carry, such as `plain`, `snv`, `mnv`, `lowqual`, `insertion`, and `deletion`.

## Development and Testing

See the [contributing guide](https://github.com/clintval/streampile/blob/main/CONTRIBUTING.md) for more information.

The streaming design follows the `StreamingPileupBuilder` of [fgbio](https://github.com/fulcrumgenomics/fgbio) of which I was also the author.

See [NOTICE](https://github.com/clintval/streampile/blob/main/NOTICE) for official attribution.
