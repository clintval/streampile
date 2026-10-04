import argparse
import sys
from pathlib import Path

from bedspec import Bed3N
from bedspec import BedReader
from bedspec import Territory
from pybgzf import IndexFormat
from pysam import AlignmentFile
from pysam import FastaFile

from streampile._table import BGZF_SUFFIXES
from streampile._table import TabulationWriter
from streampile._tabulate import DEFAULT_EXCLUDE_FLAGS
from streampile._tabulate import tabulate


def _tabulate(args: argparse.Namespace) -> int:
    with BedReader.from_path[Bed3N](args.intervals) as features:
        territory = Territory(features)
    with (
        AlignmentFile(str(args.bam), threads=args.threads) as alignments,
        FastaFile(str(args.ref)) as reference,
    ):
        bases = tabulate(
            alignments,
            reference,
            territory,
            min_base_quality=args.min_base_quality,
            min_mapping_quality=args.min_mapping_quality,
            exclude_flags=args.exclude_flags,
        )
        bgzf = args.out.suffix in BGZF_SUFFIXES
        with TabulationWriter.from_path(
            args.out,
            index=None if args.index is None else IndexFormat[args.index.upper()],
            threads=args.threads if bgzf else 1,
        ) as writer:
            writer.write_header()
            for base in bases:
                writer.write(base)
    return 0


def _flags(text: str) -> int:
    return int(text, 0)


def main(argv: list[str] | None = None) -> int:
    """Run the `streampile` command line."""
    parser = argparse.ArgumentParser(prog="streampile", description="Streaming pileups of reads.")
    commands = parser.add_subparsers(dest="command", required=True)
    tabulate_parser = commands.add_parser(
        "tabulate",
        help="count the reads of every allele at every base of a territory",
        description=(
            "Write one row for every base of a BED territory, covered or not: the reference base, "
            "the informative depth, the reference reads, and the reads of each normalized VCF "
            "allele anchored at the base, as parallel comma-separated lists."
        ),
    )
    tabulate_parser.add_argument(
        "--bam", type=Path, required=True, help="indexed, coordinate-sorted SAM, BAM, or CRAM"
    )
    tabulate_parser.add_argument("--ref", type=Path, required=True, help="indexed reference FASTA")
    tabulate_parser.add_argument(
        "--intervals", type=Path, required=True, help="the territory, a BED file"
    )
    tabulate_parser.add_argument(
        "--min-base-quality",
        type=int,
        default=0,
        help="lowest base quality of an informative base (default: 0)",
    )
    tabulate_parser.add_argument(
        "--min-mapping-quality",
        type=int,
        default=0,
        help="lowest mapping quality of a counted read (default: 0)",
    )
    tabulate_parser.add_argument(
        "--exclude-flags",
        type=_flags,
        default=DEFAULT_EXCLUDE_FLAGS,
        help=f"leave out reads with any of these SAM flags (default: {DEFAULT_EXCLUDE_FLAGS:#x})",
    )
    tabulate_parser.add_argument(
        "--index",
        choices=["tbi", "csi"],
        help="index a BGZF table as it is written (default: no index)",
    )
    tabulate_parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="threads decompressing the BAM and compressing a BGZF table (default: 1)",
    )
    tabulate_parser.add_argument(
        "--out",
        type=Path,
        required=True,
        help="the table to write, as BGZF if it ends in .gz or .bgz",
    )
    tabulate_parser.set_defaults(run=_tabulate)
    args = parser.parse_args(argv)
    return int(args.run(args))


if __name__ == "__main__":
    sys.exit(main())
