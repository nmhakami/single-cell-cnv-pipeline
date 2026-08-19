import argparse

import numpy as np
import pandas as pd
from statsmodels.nonparametric.smoothers_lowess import lowess


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate a normalized single-cell CNV profile."
    )

    parser.add_argument(
        "--counts",
        required=True,
        help="BED-like file containing genomic bins and read counts.",
    )

    parser.add_argument(
        "--gc",
        required=True,
        help="BED-like file containing GC content for each bin.",
    )

    parser.add_argument(
        "--bad-bins",
        required=False,
        default=None,
        help="Optional BED file containing bins to mask.",
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output file for the normalized CNV profile.",
    )

    parser.add_argument(
        "--lowess-fraction",
        type=float,
        default=0.3,
        help="LOWESS smoothing fraction.",
    )

    return parser.parse_args()


def read_bin_counts(path):
    """Read genomic bins and raw read counts."""
    return pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "count"],
    )


def read_gc_content(path):
    """Read GC content for each genomic bin."""
    return pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "gc"],
    )


def read_bad_bins(path):
    """
    Read genomic bins that should be masked.

    Returns an empty set when no bad-bin file is provided.
    """
    if path is None:
        return set()

    bad_bins = set()

    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            fields = line.rstrip().split("\t")

            if len(fields) < 3:
                continue

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])

            bad_bins.add((chrom, start, end))

    return bad_bins


def merge_counts_and_gc(counts, gc):
    """Match GC values to the corresponding genomic bins."""
    return counts.merge(
        gc,
        on=["chrom", "start", "end"],
        how="inner",
        validate="one_to_one",
    )


def mask_bad_bins(df, bad_bins):
    """
    Mask bins identified as problematic.

    Masked bins have their count set to NaN and are excluded
    from downstream normalization and GC correction.
    """
    if not bad_bins:
        df["bad_bin"] = False
        return df

    df = df.copy()

    df["bad_bin"] = [
        (chrom, start, end) in bad_bins
        for chrom, start, end
        in zip(df["chrom"], df["start"], df["end"])
    ]

    df.loc[df["bad_bin"], "count"] = np.nan

    return df


def normalize_read_depth(df):
    """
    Normalize bin counts for differences in sequencing depth.

    Uses the median count among unmasked bins.
    """
    valid_counts = df["count"].dropna()

    if valid_counts.empty:
        raise ValueError("No valid bins remain after masking.")

    median_count = valid_counts.median()

    if median_count <= 0:
        raise ValueError("Median bin count must be greater than zero.")

    df = df.copy()
    df["normalized_count"] = df["count"] / median_count

    return df


def gc_correct(df, fraction=0.3):
    """
    Correct GC-associated bias using LOWESS.

    The LOWESS curve models normalized read depth as a function
    of GC content among unmasked bins.
    """
    df = df.copy()

    valid = (
        ~df["bad_bin"]
        & df["gc"].notna()
        & np.isfinite(df["gc"])
        & df["normalized_count"].notna()
        & np.isfinite(df["normalized_count"])
    )

    if valid.sum() < 3:
        raise ValueError(
            "Not enough valid bins for LOWESS correction."
        )

    fitted = lowess(
        endog=df.loc[valid, "normalized_count"],
        exog=df.loc[valid, "gc"],
        frac=fraction,
        return_sorted=False,
    )

    df["gc_expected"] = np.nan
    df.loc[valid, "gc_expected"] = fitted

    df["gc_corrected"] = np.nan
    df.loc[valid, "gc_corrected"] = (
        df.loc[valid, "normalized_count"]
        / df.loc[valid, "gc_expected"]
    )

    return df


def write_profile(df, output):
    """Write the corrected bin-level CNV profile."""
    columns = [
        "chrom",
        "start",
        "end",
        "count",
        "gc",
        "bad_bin",
        "normalized_count",
        "gc_expected",
        "gc_corrected",
    ]

    df[columns].to_csv(
        output,
        sep="\t",
        index=False,
    )


def main():
    args = parse_args()

    counts = read_bin_counts(args.counts)
    gc = read_gc_content(args.gc)
    bad_bins = read_bad_bins(args.bad_bins)

    profile = merge_counts_and_gc(
        counts,
        gc,
    )

    profile = mask_bad_bins(
        profile,
        bad_bins,
    )

    profile = normalize_read_depth(
        profile,
    )

    profile = gc_correct(
        profile,
        fraction=args.lowess_fraction,
    )

    write_profile(
        profile,
        args.output,
    )


if __name__ == "__main__":
    main()