import argparse
from bisect import bisect_right

import pysam


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate genomic bin counts from a filtered BAM file."
    )

    parser.add_argument(
        "--bam",
        required=True,
        help="Path to the filtered BAM file.",
    )
    parser.add_argument(
        "--bins",
        required=True,
        help="BED file containing genomic bin coordinates.",
    )
    parser.add_argument(
        "--exclude",
        required=True,
        help="BED file containing regions to exclude.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output file for bin counts.",
    )

    return parser.parse_args()


def read_bins(bins_path):
    """
    Read genomic bins and prepare chromosome-specific lookup structures.

    Returns
    -------
    bins : list[tuple]
        Original bin definitions as (chrom, start, end).
    bins_by_chrom : dict
        Chromosome-indexed bin lookup information.
    """
    bins = []
    bins_by_chrom = {}

    with open(bins_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            fields = line.rstrip().split("\t")
            if len(fields) < 3:
                continue

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])

            bin_index = len(bins)
            bins.append((chrom, start, end))

            if chrom not in bins_by_chrom:
                bins_by_chrom[chrom] = {
                    "starts": [],
                    "ends": [],
                    "indices": [],
                }

            bins_by_chrom[chrom]["starts"].append(start)
            bins_by_chrom[chrom]["ends"].append(end)
            bins_by_chrom[chrom]["indices"].append(bin_index)

    return bins, bins_by_chrom


def read_excluded_regions(exclude_path):
    """
    Read excluded genomic regions and organize them by chromosome.
    """
    excluded_by_chrom = {}

    with open(exclude_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            fields = line.rstrip().split("\t")
            if len(fields) < 3:
                continue

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])

            excluded_by_chrom.setdefault(chrom, []).append((start, end))

    # Sort regions to allow efficient interval scanning.
    for chrom in excluded_by_chrom:
        excluded_by_chrom[chrom].sort()

    return excluded_by_chrom


def overlaps_excluded_region(position, regions):
    """
    Check whether a genomic position falls within an excluded region.

    Regions are assumed to be sorted by start coordinate.
    """
    if not regions:
        return False

    starts = [region[0] for region in regions]
    index = bisect_right(starts, position) - 1

    if index >= 0:
        start, end = regions[index]
        if start <= position < end:
            return True

    return False


def find_bin(position, bin_info):
    """
    Find the bin containing a genomic position.

    Uses binary search rather than scanning every bin.
    """
    starts = bin_info["starts"]
    ends = bin_info["ends"]
    indices = bin_info["indices"]

    index = bisect_right(starts, position) - 1

    if index >= 0 and position < ends[index]:
        return indices[index]

    return None


def count_reads(bam_path, bins, bins_by_chrom, excluded_by_chrom):
    """
    Count eligible alignments in each genomic bin.

    Reads that are unmapped, outside defined bins, or inside excluded
    regions are ignored.
    """
    counts = [0] * len(bins)

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue

            chrom = read.reference_name
            position = read.reference_start

            if chrom not in bins_by_chrom:
                continue

            if overlaps_excluded_region(
                position,
                excluded_by_chrom.get(chrom, []),
            ):
                continue

            bin_index = find_bin(
                position,
                bins_by_chrom[chrom],
            )

            if bin_index is not None:
                counts[bin_index] += 1

    return counts


def write_counts(output_path, bins, counts):
    """
    Write genomic bin coordinates and read counts to a BED-like file.
    """
    with open(output_path, "w") as f:
        for i, (chrom, start, end) in enumerate(bins):
            f.write(
                f"{chrom}\t{start}\t{end}\t{counts[i]}\n"
            )


def main():
    args = parse_args()

    bins, bins_by_chrom = read_bins(args.bins)
    excluded_by_chrom = read_excluded_regions(args.exclude)

    counts = count_reads(
        args.bam,
        bins,
        bins_by_chrom,
        excluded_by_chrom,
    )

    write_counts(
        args.output,
        bins,
        counts,
    )


if __name__ == "__main__":
    main()