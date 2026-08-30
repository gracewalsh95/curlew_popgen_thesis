#!/usr/bin/env python3
"""
Extract genes containing or within 5 kb of pcadapt outlier SNP positions.

Usage:
    python extract_genes_pcadapt.py -g annotation.gff3 -s outlier_positions.bim -o output_prefix

Input:
    -g: GFF3 annotation file
    -s: Outlier positions in .bim format (whitespace-delimited: chrom, SNP_ID, genetic_dist, position, A1, A2)

Output:
    {prefix}_outliers_annotations_autosomes.txt - All overlapping features
    {prefix}_outliers_with_genes_autosomes.txt  - Features with gene names
    {prefix}_gene_ID_for_GO_autosomes.txt       - Gene-type features only
"""

import pandas as pd
import numpy as np
import allel
import argparse

parser = argparse.ArgumentParser(
    description="Find overlaps between GFF3 and pcadapt outlier positions"
)
parser.add_argument("-g", "--gff", type=str, required=True, help="GFF3 annotation file")
parser.add_argument("-s", "--snps", type=str, required=True, help="Outlier positions (.bim format)")
parser.add_argument("-o", "--output", type=str, required=True, help="Output file prefix")
opts = parser.parse_args()

# Read outlier positions (extracted from .bim file)
outliers = pd.read_csv(
    opts.snps, sep=r"\s+", header=None,
    names=["chromosome", "SNP_ID", "genetic_dist", "position", "A1", "A2"]
)
outliers["seqid"] = outliers["chromosome"].astype(str)

annotations = allel.gff3_to_dataframe(
    opts.gff, attributes=["type", "start", "end", "ID", "Name", "biotype", "description"]
)

# Build lookup lists
start_snp = outliers["position"].tolist()
seqid_snp = outliers["seqid"].tolist()
chrom_snp = outliers["chromosome"].tolist()

start_anno = annotations.iloc[:, 3].tolist()
end_anno = annotations.iloc[:, 4].tolist()
ID_anno = annotations.iloc[:, 11].tolist()
Name_anno = annotations["Name"].tolist()
gene_anno = annotations.iloc[:, -1].tolist()
seqid_anno = annotations.iloc[:, 0].tolist()
type_anno = annotations.iloc[:, 2].tolist()

# Find annotations overlapping each outlier SNP position
BUFFER = 5000  # ±5 kb of gene body

results = []
for n in range(len(start_snp)):
    for i in range(len(seqid_anno)):
        if (seqid_snp[n] == seqid_anno[i]) and \
           (start_anno[i] >= (start_snp[n] - BUFFER)) and \
           (end_anno[i] <= (start_snp[n] + BUFFER)):
            results.append([
                chrom_snp[n], str(start_snp[n]), str(start_snp[n]),
                seqid_anno[i], str(start_anno[i]), str(end_anno[i]),
                type_anno[i], ID_anno[i], Name_anno[i], gene_anno[i]
            ])

annotation_and_snps = pd.DataFrame(results, columns=[
    "chrom", "start_bin", "end_bin", "seqID",
    "start_anno", "end_anno", "type_anno", "ID", "Name", "gene"
])

annotation_and_snps.to_csv(opts.output + "_outliers_annotations_autosomes.txt", sep="\t")

genes = annotation_and_snps[annotation_and_snps.iloc[:, -1] != "."]
genes.to_csv(opts.output + "_outliers_with_genes_autosomes.txt", sep="\t")

ID_for_GO = annotation_and_snps[annotation_and_snps["type_anno"] == "gene"]
ID_for_GO.to_csv(opts.output + "_gene_ID_for_GO_autosomes.txt", sep="\t")

print(f"Outlier positions: {len(outliers)}")
print(f"Gene annotations found: {len(ID_for_GO)}")
print(f"Unique genes: {ID_for_GO['Name'].nunique()}")