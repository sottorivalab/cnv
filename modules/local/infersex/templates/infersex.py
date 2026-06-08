#!/usr/bin/env python

from pathlib import Path
from statistics import median

import pysam


NORMAL_BAM = "${normal_bam}"
PREFIX = "${prefix}"

REGIONS = {
    "autosomes": [
        ("chr1", 10_000_000, 10_100_000),
        ("chr2", 20_000_000, 20_100_000),
        ("chr3", 30_000_000, 30_100_000),
        ("chr4", 40_000_000, 40_100_000),
        ("chr5", 50_000_000, 50_100_000),
        ("chr6", 60_000_000, 60_100_000),
    ],
    "chrX": [
        ("chrX", 50_000_000, 50_100_000),
        ("chrX", 60_000_000, 60_100_000),
        ("chrX", 70_000_000, 70_100_000),
        ("chrX", 80_000_000, 80_100_000),
        ("chrX", 90_000_000, 90_100_000),
        ("chrX", 100_000_000, 100_100_000),
    ],
    "chrY": [
        ("chrY", 7_000_000, 7_100_000),
        ("chrY", 8_000_000, 8_100_000),
        ("chrY", 9_000_000, 9_100_000),
        ("chrY", 13_000_000, 13_100_000),
        ("chrY", 14_000_000, 14_100_000),
        ("chrY", 15_000_000, 15_100_000),
    ],
}


def resolve_chromosome_name(bam, chrom):
    references = set(bam.references)
    if chrom in references:
        return chrom
    without_chr = chrom.removeprefix("chr")
    if without_chr in references:
        return without_chr
    raise ValueError(f"Chromosome {chrom!r} was not found in the BAM header")


def mean_depth(bam, chrom, start, end):
    resolved_chrom = resolve_chromosome_name(bam, chrom)
    coverage = bam.count_coverage(resolved_chrom, start, end)
    depths = [sum(base_depths) for base_depths in zip(*coverage)]
    return sum(depths) / len(depths)


def identify_sex(x_depth, y_depth, autosome_depth):
    if autosome_depth <= 0:
        return "unknown"

    x_ratio = x_depth / autosome_depth
    y_ratio = y_depth / autosome_depth

    if x_ratio <= 0.65 and y_ratio >= 0.10:
        return "male"
    if x_ratio >= 0.65 and y_ratio <= 0.10:
        return "female"
    return "unknown"


rows = []

with pysam.AlignmentFile(NORMAL_BAM, "rb") as bam:
    for group, regions in REGIONS.items():
        for chrom, start, end in regions:
            depth = mean_depth(bam, chrom, start, end)
            rows.append((group, chrom, start, end, depth))

group_depths = {}
for group in REGIONS:
    group_depths[group] = median(row[4] for row in rows if row[0] == group)

autosome_depth = group_depths["autosomes"]
x_depth = group_depths["chrX"]
y_depth = group_depths["chrY"]
sex = identify_sex(x_depth, y_depth, autosome_depth)

Path(f"{PREFIX}.sex.txt").write_text(f"{sex}\n")

with Path(f"{PREFIX}.sex_depths.tsv").open("w") as handle:
    handle.write("group\tchrom\tstart\tend\tdepth\n")
    for group, chrom, start, end, depth in rows:
        handle.write(f"{group}\t{chrom}\t{start}\t{end}\t{depth:.6f}\n")

print(f"autosome median: {autosome_depth:.6f}")
print(f"chrX/autosome: {x_depth / autosome_depth if autosome_depth else 'NA'}")
print(f"chrY/autosome: {y_depth / autosome_depth if autosome_depth else 'NA'}")
print(f"inferred sex: {sex}")

with Path("versions.yml").open("w") as handle:
    handle.write('"${task.process}":\n')
    handle.write(f"    pysam: {pysam.__version__}\n")
