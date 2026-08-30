# Chapter 4: Inbreeding, ROH, and Genetic Load

Code and parameters to reproduce the results in Chapter 4. Standard bioinformatics pipeline steps (read trimming, mapping, variant calling, quality filtering) are described in Sections 4.2.1 and Chapter 3 Sections 3.2.1–3.2.4 and follow GATK best practices. The input VCF (no-MAF filtered, ~7.5M autosomal biallelic SNPs) was generated in Chapter 3.

**Reference genome:** GCA_964106895.1 (*Numenius arquata*)

**Samples:** 56 after QC (Ireland 23, Britain 18, Sweden 15). Curl1 (reference genome individual) excluded from per-individual estimates.

## Scripts

| Script | Description |
|--------|-------------|
| `build_load_table.sh` | Splits SnpEff-annotated VCF by impact class (LOW/MODERATE/HIGH), builds per-individual table of masked (het) and realised (hom) derived deleterious genotype counts, appends population labels |
| `rxy_jackknife.R` | Calculates Rxy (Do et al. 2015) for each population pair, standardised to non-coding (MODIFIER) sites. CIs via weighted block jackknife (Busing et al. 1999) over ~115 contiguous blocks. Block generation and Busing SE functions from [Stuart (2026)](https://github.com/OliverPStuart/2025_Jackknife_Review); SNP-to-block assignment vectorised with `data.table::foverlaps` (identical assignments) |
| `roh_analysis.R` | Processes BCFtools/RoH output into short (100 kb – 1 Mb) and long (> 1 Mb) ROH. Calculates F_ROH. Kruskal-Wallis and pairwise Wilcoxon tests (BH correction). Spearman correlations of per-window (1 Mb) ROH proportions. Generates heatmaps, violin plots, scatter plots |
| `outlier_roh_load_overlap.sh` | Three-way intersection of Ch3 selection outlier genes, long ROH, and SnpEff impact annotations. Maps exon-level transcript IDs to parent genes. Extracts per-individual genotypes for classifying homozygous derived genotypes within/outside ROH |
| `plot_venn.R` | Venn diagrams of shared/private deleterious variants from `bcftools isec` output. Decodes positional binary CODE column (Brit=pos1, Ire=pos2, Swe=pos3) |

## Analysis Parameters

### HWE filtering

VCFtools v0.1.16. SNPs deviating at *P* < 1 × 10⁻⁴ in any population excluded for ROH analysis.

```bash
vcftools --gzvcf input.vcf.gz --hwe 0.0001 --recode --out output
```

### Allele frequency estimation for ROH

BCFtools v1.15.1. Allele frequencies estimated from study samples, formatted for BCFtools/RoH input.

```bash
bcftools +fill-tags filtered.vcf.gz -- -t AF | \
  bcftools query -f '%CHROM\t%POS\t%REF,%ALT\t%AF\n' | \
  sort -k1,1 -k2,2n > af.txt

sed -i '1i #CHROM\tPOS\tREF,ALT\tAF' af.txt
bgzip -f af.txt
tabix -f -s1 -b2 -e2 af.txt.gz
```

### ROH calling

BCFtools/RoH v1.15.1 (Narasimhan et al. 2016). HMM on genotype likelihoods. Recombination rate 2.8 cM/Mb for avian genomes. ROH defined as tracts > 100 kb, split into short (100 kb – 1 Mb) and long (> 1 Mb) in `roh_analysis.R`.

```bash
bcftools roh \
  --rec-rate 2e-8 \
  --AF-file af.txt.gz \
  --output-type r \
  -o output.txt \
  filtered_hwe.vcf.gz
```

### Ancestral allele polarisation

Outgroup species for ancestral state reconstruction:

| Species | Accession | Data |
|---------|-----------|------|
| Little curlew (*N. minutus*) | SRX12367832 | 30× WGS |
| Whimbrel (*N. phaeopus*) | GCA_030770645.1 | Scaffold-level assembly |
| Upland sandpiper (*B. longicauda*) | ASM4578454v1 | Scaffold-level assembly |

Reference assemblies processed following Kardos (2023): masked sequences removed, N positions as breakpoints, remaining sequences split into uniform lengths, simulated as 100 bp overlapping reads at 10 bp intervals using BEDTools getfasta, mapped to curlew reference with BWA-mem. Little curlew WGS mapped directly. See [Kardos (2023)](https://github.com/martinkardos/WolfGenomePurging).

**Consensus calling** with ANGSD v0.940:

```bash
angsd -bam bam_list.txt \
  -doFasta 2 -doCounts 1 -P 8 \
  -rf autosomes_list.txt \
  -out outgroup_consensus
```

**Consensus parsing** (extract polarisable positions and alleles):

```bash
awk '
  NR==FNR {autosomes[$1]=1; next}
  /^>/ {chrom=substr($0,2); pos=0; keep=(chrom in autosomes); next}
  keep {
    n=split($0, bases, "")
    for(i=1; i<=n; i++) {
      pos++; base=toupper(bases[i])
      if(base != "N") {
        print chrom"\t"pos >> "consensus_positions.txt"
        print chrom":"pos"\t"base >> "consensus_alleles.txt"
      }
    }
  }
' autosomes.txt outgroup_consensus.fa
```

**Chromosome name conversions:** Two mapping files bridge three naming schemes: numeric (PLINK-friendly) → ENA accession (consensus/polarisation) → NC_ (SnpEff database).

**Step 1:** Rename chromosomes (numeric → ENA):

```bash
bcftools annotate --rename-chrs correct_chrom_map.txt input.vcf.gz -Oz -o ena_names.vcf.gz
```

**Step 2:** Extract consensus SNP positions with PLINK2 v2.0.0:

```bash
awk '{print $1":"$2}' consensus_positions.txt > consensus_positions_fixed.txt

plink2 --vcf ena_names.vcf.gz \
  --allow-extra-chr --chr-set 40 \
  --set-all-var-ids @:# \
  --extract consensus_positions_fixed.txt \
  --make-pgen --out curlew_filtered
```

**Step 3:** Polarise (force REF = ancestral allele):

```bash
sed 's/\t/ /g' consensus_alleles.txt > consensus_alleles_fixed.txt

plink2 --pfile curlew_filtered \
  --allow-extra-chr --chr-set 40 \
  --ref-allele force consensus_alleles_fixed.txt \
  --make-pgen --out curlew_polarized
```

**Step 4:** Extract mismatches from polarisation log (must be regenerated per run):

```bash
grep 'Warning' curlew_polarized.log | \
  awk '{print $7}' | sed "s/'//g" | sed "s/\.$//g" > mismatches.txt
```

**Step 5:** Export to VCF excluding mismatches:

```bash
plink2 --pfile curlew_polarized \
  --exclude mismatches.txt \
  --allow-extra-chr --chr-set 40 \
  --export vcf-4.2 bgz \
  --out polarized_final
```

**Step 6:** Rename ENA → NC_ for SnpEff database:

```bash
bcftools annotate --rename-chrs final_ena_nc_bridge.txt \
  polarized_final.vcf.gz -Oz -o polarized_nc.vcf.gz
```

### SnpEff annotation

SnpEff v5.2c. Custom database (bNumArq3) built from GCA_964106895.1 annotation. Canonical transcripts only.

```bash
java -Xmx8g -jar snpEff.jar ann -v -canon \
  -c snpEff.config \
  -stats report.html \
  bNumArq3 \
  polarized_nc.vcf.gz 2> snpeff.log | bgzip > annotated_autosomes.vcf.gz
```

### Rxy data preparation

Site classification and per-population allele frequencies extracted from annotated VCF for input to `rxy_jackknife.R`.

```bash
SnpSift extractFields annotated_autosomes.vcf.gz \
  CHROM POS "ANN[0].IMPACT" "ANN[0].EFFECT" > site_class.txt

for pop in Ire_nocurl1 Brit Swe; do
  bcftools view -S ${pop}.txt --force-samples annotated_autosomes.vcf.gz | \
    bcftools +fill-tags -- -t AF | \
    bcftools query -f '%CHROM\t%POS\t%INFO/AF\n' > af_${pop}.txt
done
```

### Degenerate site diversity (π₀/π₄)

Degeneracy annotation with degenotate v1.0 (Mirchandani et al. 2024). 0-fold and 4-fold sites extracted to BED files, intersected with per-chromosome all-sites VCFs.

```bash
degenotate -g annotation.gff -f reference.fa -o degeneracy_scores.txt

bedtools intersect \
  -a chr_allsites.vcf.gz \
  -b 0fold_sorted.bed.gz \
  -header | bgzip > chr_0fold.vcf.gz
```

Nucleotide diversity with pixy v1.2.7 (Korunes & Samuk 2021):

```bash
pixy --stats pi \
  --vcf chr_0fold.vcf.gz \
  --populations all_3popfile.tsv \
  --window_size 15000 \
  --n_cores 6 \
  --output_folder output
```

Run per chromosome (1–39), for both 0-fold and 4-fold, at window sizes 15,000 bp and 25,000 bp. The ratio π₀/π₄ was calculated per population as an estimate of purifying selection efficacy.

### Venn diagram (shared/private deleterious variants)

Per-population variant extraction and intersection with BCFtools v1.15.1 and SnpSift. Population order matters for CODE decoding.

```bash
for POP in Brit Ire Swe; do
  bcftools view -S ${POP}.txt --min-ac 1 annotated_autosomes.vcf.gz \
    -Oz -o ${POP}_annotated.vcf.gz

  SnpSift filter "(ANN[0].IMPACT = 'HIGH') | (ANN[0].IMPACT = 'MODERATE')" \
    ${POP}_annotated.vcf.gz | bgzip > ${POP}_high_mod.vcf.gz

  SnpSift filter "(ANN[0].IMPACT = 'HIGH')" \
    ${POP}_annotated.vcf.gz | bgzip > ${POP}_high.vcf.gz

  SnpSift filter "(ANN[0].IMPACT = 'MODERATE')" \
    ${POP}_annotated.vcf.gz | bgzip > ${POP}_mod.vcf.gz
done

# Intersection (order: Brit, Ire, Swe → CODE positions 1, 2, 3)
bcftools isec -n +1 -p intersection \
  Brit_high_mod.vcf.gz Ire_high_mod.vcf.gz Swe_high_mod.vcf.gz

bcftools isec -n +1 -p intersection/high \
  Brit_high.vcf.gz Ire_high.vcf.gz Swe_high.vcf.gz

bcftools isec -n +1 -p intersection/mod \
  Brit_mod.vcf.gz Ire_mod.vcf.gz Swe_mod.vcf.gz
```

## Notes

Chromosome name mapping files (`correct_chrom_map.txt`, `final_ena_nc_bridge.txt`) bridge three naming schemes: numeric (PLINK-friendly), ENA accession (consensus/polarisation), and NC_ (SnpEff database).

