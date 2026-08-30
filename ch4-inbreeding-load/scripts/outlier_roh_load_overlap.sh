#!/bin/bash
# outlier_roh_load_overlap.sh — intersect Chapter 3 outlier genes with long ROH
#                                and SnpEff impact annotations
#
# Usage:
# bash outlier_roh_load_overlap.sh <outlier_list1> <outlier_list2> \
#      <long_roh_file> <annotated_vcf> <snpeff_genes_txt> \
#      <exclude_sample> <output_dir>
#
# Arguments:
#   outlier_list1      TSV with columns: chrom, region, region_start, region_end,
#                      start_anno, end_anno, ID, Name, gene (header line)
#   outlier_list2      TSV with columns: chrom, start_bin, end_bin,
#                      start_anno, end_anno, ID, Name, gene (header line)
#   long_roh_file      Long ROH file (>1 Mb), with header. Columns: RG, Sample,
#                      Chromosome, Start, End, Length_bp, Num_markers, Quality,
#                      Pop, KB, ROH_category
#   annotated_vcf      SnpEff-annotated VCF (bgzipped + indexed)
#   snpeff_genes_txt   SnpEff gene summary file (e.g. Curlew_Final_Report.genes.txt)
#   exclude_sample     Sample ID to exclude (e.g. Curl1)
#   output_dir         Output directory
#
# Requires: bedtools, bcftools, tabix, bgzip, awk, grep, sed
#
# Output:
#   outlier_in_roh_with_moderate.txt   Genes at intersection of all three lists for predicted moderate impact variants
#   outlier_in_roh_with_high.txt       Genes at intersection of all three lists for predicted high impact variants)
#   genotypes_only.txt                 Per-individual genotypes for R import
#   sample_names.txt                   Sample order
#   variant_genes_mapped.txt           Variants mapped to gene names


OUTLIER1="$1"
OUTLIER2="$2"
LONG_ROH="$3"
ANNO_VCF="$4"
GENES_TXT="$5"
EXCLUDE="$6"
OUT="$7"

mkdir -p "$OUT"
cd "$OUT"

# 1. Convert outlier gene lists to BED (I used two methods to detect outliers so had two lists)
tail -n +2 "$OUTLIER1" | awk -F'\t' '{print $1, $5-1, $6, $7, $9}' OFS='\t' > outlier1.bed
tail -n +2 "$OUTLIER2" | awk -F'\t' '{print $1, $4-1, $5, $6, $8}' OFS='\t' > outlier2.bed
cat outlier1.bed outlier2.bed | sort -k1,1 -k2,2n | uniq > all_outlier_genes.bed

# 2. Convert long ROH to BED
tail -n +2 "$LONG_ROH" | awk -F'\t' '{print $3, $4-1, $5, $2, $9}' OFS='\t' > long_roh.bed

# 3. Intersect outlier genes with long ROH
bedtools intersect -a all_outlier_genes.bed -b long_roh.bed -wa -wb > outlier_genes_in_long_roh.bed

# 4. Extract unique gene names from intersection
awk -F'\t' '$5 != "." {print $5}' outlier_genes_in_long_roh.bed | sort -u > outlier_roh_gene_names.txt

# 5. Find protein-coding genes with moderate-impact variants
# Map exon transcript IDs back to parent genes via the SnpEff gene summary
awk -F'\t' '$7 > 0 && $1 ~ /^exon-XM/' "$GENES_TXT" | \
  sed 's/exon-\(XM_[0-9]*\.[0-9]*\).*/\1/' | \
  sort -u > transcripts_with_moderate.txt

awk -F'\t' '$4 == "protein_coding" {print $3, $1}' OFS='\t' "$GENES_TXT" | \
  grep -Fwf transcripts_with_moderate.txt | \
  awk -F'\t' '{print $2}' | sort -u > moderate_protein_coding_genes.txt

# 6. Same for high-impact variants (column 5)
awk -F'\t' '$5 > 0 && $1 ~ /^exon-XM/' "$GENES_TXT" | \
  sed 's/exon-\(XM_[0-9]*\.[0-9]*\).*/\1/' | \
  sort -u > transcripts_with_high.txt

awk -F'\t' '$4 == "protein_coding" {print $3, $1}' OFS='\t' "$GENES_TXT" | \
  grep -Fwf transcripts_with_high.txt | \
  awk -F'\t' '{print $2}' | sort -u > high_protein_coding_genes.txt

# 7. Three-way intersection: outlier + ROH + impact
grep -iFxf outlier_roh_gene_names.txt moderate_protein_coding_genes.txt \
  > outlier_in_roh_with_moderate.txt || true
grep -iFxf outlier_roh_gene_names.txt high_protein_coding_genes.txt \
  > outlier_in_roh_with_high.txt || true

echo "Moderate-impact outlier genes in long ROH: $(wc -l < outlier_in_roh_with_moderate.txt)"
echo "High-impact outlier genes in long ROH: $(wc -l < outlier_in_roh_with_high.txt)"

# 8. Extract moderate-impact variants for outlier genes
if [ -s outlier_in_roh_with_moderate.txt ]; then
  zgrep "^#" "$ANNO_VCF" > full_header.txt
  zgrep -wf outlier_in_roh_with_moderate.txt "$ANNO_VCF" | grep "MODERATE" > moderate_outlier_variants.vcf
  cat full_header.txt moderate_outlier_variants.vcf > moderate_outlier_with_header.vcf
  bgzip -f moderate_outlier_with_header.vcf
  tabix -p vcf moderate_outlier_with_header.vcf.gz

  # Exclude reference-genome individual
  bcftools view -s "^${EXCLUDE}" moderate_outlier_with_header.vcf.gz \
    -Oz -o moderate_outlier_filtered.vcf.gz
  tabix -p vcf moderate_outlier_filtered.vcf.gz

  #9. Extract genotypes for R import
  bcftools query -f '%CHROM\t%POS[\t%GT]\n' moderate_outlier_filtered.vcf.gz > genotypes_only.txt
  bcftools query -l moderate_outlier_filtered.vcf.gz > sample_names.txt

  # 10. Map variants to gene names via transcript IDs
  bcftools query -f '%CHROM\t%POS\t%INFO/ANN\n' moderate_outlier_filtered.vcf.gz | \
    awk -F'\t' '{split($3,a,"|"); print $1"\t"$2"\t"a[2]"\t"a[4]}' > variant_info.txt

  awk -F'\t' '{match($4, /XM_[0-9]+\.[0-9]+/, m); print $0"\t"m[0]}' \
    variant_info.txt > variant_with_transcript.txt

  awk -F'\t' '$4 == "protein_coding" {print $3, $1}' OFS='\t' \
    "$GENES_TXT" > transcript_to_gene.txt

  awk -F'\t' 'NR==FNR {gene[$1]=$2; next} {print $0"\t"gene[$5]}' \
    transcript_to_gene.txt variant_with_transcript.txt > variant_genes_mapped.txt
fi

echo "Done. Output in ${OUT}/"
