#!/bin/bash
#
# Usage: bash build_load_table.sh <annotated_vcf> <population_file> <output_dir>
#
# Arguments:
# annotated_vcf    SnpEff-annotated VCF (bgzipped), REF = ancestral allele
# population_file  Tab-delimited file: sample_id <tab> population_name
# output_dir       Directory for output files
#
#
# Output:
#   low_impact.vcf
#   moderate_impact.vcf
#   high_impact.vcf
#   Load_table.txt              (all individuals)
#   Load_table_with_pop.txt      (with population column)
#
# Load definitions:
#   masked   = heterozygous derived deleterious genotypes (0/1 or 0|1)
#   realised = homozygous derived deleterious genotypes (1/1 or 1|1)

ANNOTATED_VCF="$1"
POP_FILE="$2"
OUT="$3"

mkdir -p "$OUT"
cd "$OUT"

# Split by impact class
for impact in LOW MODERATE HIGH; do
  out_vcf=$(echo "$imp" | tr '[:upper:]' '[:lower:]')_impact.vcf
  echo "Extracting ${imp} impact variants..."
  SnpSift filter "ANN[*].IMPACT = '${imp}'" "$ANNOTATED_VCF" > "$out_vcf"
done

# Build per-individual load table
echo -e "Type\tInd\tLoad\tNumber" > Load_table.txt

for type in low moderate high; do
  vcf_file="${type}_impact.vcf"
  header=$(grep "^#CHROM" "$vcf_file") || { echo "Error: header not found in ${vcf_file}"; exit 1; }

  for col in $(seq 10 $(echo "$header" | awk '{print NF}')); do
    ind=$(echo "$header" | awk -v c="$col" '{print $c}')
    grep -v "^#" "$vcf_file" | cut -f "$col" | cut -d ':' -f 1 | \
    awk -v OFS="\t" -v t="$type" -v i="$ind" \
      '{if($1=="0/1" || $1=="0|1") {masked++}
        else if($1=="1/1" || $1=="1|1") {realized++}}
      END{print t, i, "masked",   (masked   ? masked   : 0)
          print t, i, "realized", (realized ? realized : 0)}' >> Load_table.txt
  done
done

# Append population labels
awk 'NR==FNR {pop[$1]=substr($0, index($0,$2)); next}
     FNR==1  {print $0 "\tPopulation"}
     FNR>1   {print $0 "\t" pop[$2]}' \
  "$POP_FILE" Load_table.txt | sed 's/\t\t/\t/g' > Load_table_with_pop.txt

echo "Done. Output in ${OUT}/"
