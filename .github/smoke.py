"""Check that an installed streampile imports, reads pysam records directly, and piles them up."""

import importlib
import pkgutil

import pysam

import streampile
from streampile import StreamingPileupBuilder
from streampile._native import direct_bridge

for module in pkgutil.walk_packages(streampile.__path__, "streampile."):
    importlib.import_module(module.name)

header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 100}]}
with pysam.AlignmentFile("smoke.bam", "wb", header=header) as sink:
    read = pysam.AlignedSegment(sink.header)
    read.query_name = "r"
    read.reference_id = 0
    read.reference_start = 10
    read.mapping_quality = 60
    read.cigarstring = "2M1D2M"
    read.query_sequence = "ACGT"
    read.query_qualities = pysam.qualitystring_to_array("IIII")
    sink.write(read)
with pysam.AlignmentFile("smoke.bam") as reads, StreamingPileupBuilder(reads) as builder:
    depths = [pileup.unfiltered_depth for pileup in builder.columns("chr1", 9, 16)]
assert direct_bridge(), "pysam records are read through their attributes"
assert depths == [0, 1, 1, 1, 1, 1, 0], depths
print("streampile works")
