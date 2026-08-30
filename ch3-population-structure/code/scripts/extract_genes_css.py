#!/usr/bin/env python3
"""
Extract genes within or near CSS outlier regions from a GFF3 annotation.

Usage:
    python extract_genes_css.py -g annotation.gff3 -r merged_regions.txt -o output_prefix -b 200000

Input:
    -g: GFF3 annotation file
    -r: Tab-delimited regions file (columns: chromosome, start, end)
    -b: Buffer distance in bp added to each region boundary

Output:
    {prefix}_all_annotations.txt  - All annotation features overlapping regions
    {prefix}_genes.txt            - Features with assigned gene names
    {prefix}_genes_for_GO.txt     - Gene-type features only (for GO analysis)
"""

import pandas as pd
import numpy as np
import allel
import argparse

parser = argparse.ArgumentParser(description="Extract genes within selection regions ± buffer")
parser.add_argument("-g", "--gff", type=str, required=True, help="GFF3 annotation file")
parser.add_argument("-r", "--regions", type=str, required=True, help="Regions file (tab-delimited)")
parser.add_argument("-o", "--output", type=str, required=True, help="Output file prefix")
parser.add_argument("-b", "--buffer", type=int, required=True, help="Buffer distance in bp")
opts = parser.parse_args()

regions = pd.read_csv(opts.regions, sep="\t")
annotations = allel.gff3_to_dataframe(
    opts.gff, attributes=["type", "start", "end", "ID", "Name", "biotype", "description"]
)

regions["seqid"] = regions["chromosome"].astype(str)

# Build lookup lists
start_reg = regions["start"].tolist()
end_reg = regions["end"].tolist()
seqid_reg = regions["seqid"].tolist()
chrom_reg = regions["chromosome"].tolist()

start_anno = annotations.iloc[:, 3].tolist()
end_anno = annotations.iloc[:, 4].tolist()
ID_anno = annotations.iloc[:, 11].tolist()
Name_anno = annotations["Name"].tolist()
gene_anno = annotations.iloc[:, -1].tolist()
seqid_anno = annotations.iloc[:, 0].tolist()
type_anno = annotations.iloc[:, 2].tolist()

# Find annotations overlapping regions ± buffer
results = []
for i in range(len(seqid_anno)):
    for n in range(len(seqid_reg)):
        if (seqid_reg[n] == seqid_anno[i]) and \
           (start_anno[i] >= (start_reg[n] - opts.buffer)) and \
           (end_anno[i] <= (end_reg[n] + opts.buffer)):
            results.append([
                chrom_reg[n], str(start_reg[n]), str(end_reg[n]),
                seqid_anno[i], str(start_anno[i]), str(end_anno[i]),
                type_anno[i], ID_anno[i], Name_anno[i], gene_anno[i]
            ])

annotation_and_regions = pd.DataFrame(results, columns=[
    "chrom", "region_start", "region_end", "seqID",
    "start_anno", "end_anno", "type_anno", "ID", "Name", "gene"
])

genes = annotation_and_regions[annotation_and_regions.iloc[:, -1] != "."]
genes_only = annotation_and_regions[annotation_and_regions["type_anno"] == "gene"]

annotation_and_regions.to_csv(opts.output + "_all_annotations.txt", sep="\t", index=False)
genes.to_csv(opts.output + "_genes.txt", sep="\t", index=False)
genes_only.to_csv(opts.output + "_genes_for_GO.txt", sep="\t", index=False)

print(f"Total regions: {len(regions)}")
print(f"Gene annotations found: {len(genes_only)}")
print(f"Unique genes: {genes_only['Name'].nunique()}")
