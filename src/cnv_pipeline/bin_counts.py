import argparse
from pathlib import Path

import pysam

def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate genomic bin counts from a filtered BAM file."
    )

    parser.add_argument(
        "--bam",
        required=True,
        help="Path to the filtered BAM file."
    )

    parser.add_argument(
        "--bins",
        required=True,
        help="BED file containing genomic bin coordinates."
    )

    parser.add_argument(
        "--exclude",
        required=True,
        help="BED file containing regions to exclude."
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output file for bin counts."
    )

    return parser.parse_args()

def read_bins(bins_path):
    bins = []

    with open(bins_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            fields = line.rstrip().split("\t")

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])

            bins.append((chrom, start, end))

    return bins

def read_excluded_regions(exclude_path):
    regions = []

    with open(exclude_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            fields = line.rstrip().split("\t")

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])

            regions.append((chrom, start, end))

    return regions

def count_reads(bam_path, bins, excluded_regions):
    counts = {i: 0 for i in range(len(bins))}

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue

            chrom = read.reference_name
            start = read.reference_start

            # Skip reads that fall in an excluded region.
            excluded = any(
                region_chrom == chrom
                and region_start <= start < region_end
                for region_chrom, region_start, region_end in excluded_regions
            )

            if excluded:
                continue

            # Find the bin containing the read.
            for i, (bin_chrom, bin_start, bin_end) in enumerate(bins):
                if bin_chrom == chrom and bin_start <= start < bin_end:
                    counts[i] += 1
                    break

    return counts

def write_counts(output_path, bins, counts):
    with open(output_path, "w") as f:
        for i, (chrom, start, end) in enumerate(bins):
            f.write(f"{chrom}\t{start}\t{end}\t{counts[i]}\n")

def main():
    args = parse_args()

    bins = read_bins(args.bins)
    excluded_regions = read_excluded_regions(args.exclude)

    counts = count_reads(
        args.bam,
        bins,
        excluded_regions,
    )

    write_counts(
        args.output,
        bins,
        counts,
    )


if __name__ == "__main__":
    main()