#!/usr/bin/env python3
"""
Smooth per-SNP CSS scores over non-overlapping windows.

Usage:
    python smoothing_window.py -i CSS.txt -o CSS_smoothed_20kb.txt -w 20000

Input: tab-delimited file with columns: chromosome, position, CSS
Output: same file with added CSS_smoothed and window_id columns
"""

import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(description="Smooth CSS scores over non-overlapping windows")
parser.add_argument("-i", "--input", type=str, required=True, help="Input CSS file (tab-delimited)")
parser.add_argument("-o", "--output", type=str, required=True, help="Output file name")
parser.add_argument("-w", "--window", type=int, default=20000, help="Window size in bp (default: 20000)")
opts = parser.parse_args()

data = pd.read_csv(opts.input, sep="\t")

required = ["chromosome", "position", "CSS"]
if not all(col in data.columns for col in required):
    raise ValueError(f"Input file must contain columns: {required}")

data = data.sort_values(["chromosome", "position"]).reset_index(drop=True)
data["CSS_smoothed"] = np.nan
data["window_id"] = np.nan

for chrom in sorted(data["chromosome"].unique()):
    chr_mask = data["chromosome"] == chrom
    chr_data = data.loc[chr_mask].copy()

    # Assign each SNP to a non-overlapping bin
    chr_data["bin"] = (chr_data["position"] // opts.window).astype(int)
    bin_means = chr_data.groupby("bin")["CSS"].mean()

    chr_data["CSS_smoothed"] = chr_data["bin"].map(bin_means)
    chr_data["window_id"] = chr_data["bin"]

    data.loc[chr_mask, "CSS_smoothed"] = chr_data["CSS_smoothed"].values
    data.loc[chr_mask, "window_id"] = chr_data["window_id"].values

data.to_csv(opts.output, sep="\t", index=False)

print(f"Total SNPs: {len(data)}")
print(f"Total windows: {data['window_id'].nunique()}")
print(f"Raw CSS range: {data['CSS'].min():.3f} – {data['CSS'].max():.3f}")
print(f"Smoothed CSS range: {data['CSS_smoothed'].min():.3f} – {data['CSS_smoothed'].max():.3f}")
