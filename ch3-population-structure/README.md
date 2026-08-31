# Chapter 3: Population Structure and Adaptive Divergence

Code and parameters to reproduce the population genomics analyses in Chapter 3. Only analysis-specific scripts and non-standard parameters are included; standard pipeline steps (read trimming, mapping, GATK variant calling) follow published best practices and are described in the Methods.

**Reference genome:** GCA_964106895.1 ([doi:10.12688/wellcomeopenres.24272.1](https://doi.org/10.12688/wellcomeopenres.24272.1))

## VCFs

- **SNP VCF (MAC-filtered):** Biallelic SNPs, MAC ≥ 3, ≤ 10% missing, depth 10–40×. Used for population structure, selection scans, inbreeding, Ne.
- **All-sites VCF (no MAF filter):** Used for pixy (π, θ_W, Tajima's D, F_ST, D_xy).

## Scripts

| Script | Description |
|--------|-------------|
| `smoothing_window.py` | Averages per-SNP CSS scores in non-overlapping 20 kb windows |
| `css_outlier_detection.R` | Identifies CSS outlier SNPs (top 0.1%), validates with flanking support (≥ 5 SNPs in top 1% within ± 500 kb), defines and merges regions (≤ 200 kb) |
| `extract_genes_css.py` | Extracts genes within or ± 250 kb of CSS outlier regions from GFF3 annotation |
| `pcadapt_outliers.R` | Runs pcadapt (K = 2), applies FDR 0.01 (Storey q-value), exports outlier SNP indices |
| `extract_genes_pcadapt.py` | Extracts genes containing or ± 5 kb of pcadapt outlier SNP positions from GFF3 |
| `tajima_d_jackknife.R` | Block jackknife (1 Mb non-overlapping blocks) with pairwise t-tests for Tajima's D differences between populations |

## Analysis Parameters

### Quality control and filtering

GATK v4.3.0.0 best practices (Van Der Auwera et al. 2013). VariantFiltration: `QUAL < 30, MQ < 40, SOR > 4, QD < 2, FS > 60, MQRankSum < -12.5, ReadPosRankSum < -8`. SNPs within 5 bp of an indel removed with BEDTools v2.27.1.

```bash
vcftools --gzvcf input.vcf.gz --mac 3 --max-missing 0.9 --min-meanDP 10 --max-meanDP 40 --recode
```

### Kinship

KING v2.2.7. First-degree relatives (kinship coefficient > 0.2) removed.

```bash
king -b input.bed --kinship --related
```

### LD pruning

PLINK v1.90b6.25.

```bash
plink --vcf input.vcf --allow-extra-chr --chr-set 40 --set-missing-var-ids @:# --indep-pairwise 50 10 0.2
plink --vcf input.vcf --allow-extra-chr --chr-set 40 --set-missing-var-ids @:# --extract pruned.prune.in --make-bed --recode vcf-iid
```

### Phasing

SHAPEIT4 v4.2.2. Default recombination rate 1 cM/Mb.

```bash
shapeit4 --input input.vcf.gz --output phased.vcf.gz --region CHR
```

### PCA

PLINK v1.90b6.25. Full SNP dataset (not LD-pruned).

```bash
plink --vcf input.vcf --allow-extra-chr --chr-set 40 --set-missing-var-ids @:# --make-bed --pca
```

### ADMIXTURE

v1.3.0, LD-pruned dataset.

```bash
for k in {1..5}; do admixture --cv input_pruned.bed $k > log${k}.out; done
```

### Genetic diversity

pixy v2.0.0.beta8. Per chromosome, all-sites VCF.

```bash
pixy --stats pi watterson_theta tajima_d --vcf allsites_chrN.vcf.gz --populations popfile.tsv --window_size 15000 --n_cores 1
```

Concatenation:

```bash
cat chr1/pixy_pi.txt > combined_pixy_pi.txt
tail -n +2 chr{2..39}/pixy_pi.txt >> combined_pixy_pi.txt
```

- Genome-wide π = sum(count_diffs) / sum(count_comparisons)
- Genome-wide θ_W = weighted.mean(avg_watterson_theta, no_sites)

### F_ST and D_xy

pixy, 15 kb windows, all-sites VCF. All pairwise comparisons including merged Ireland + Britain.

### Heterozygosity

```bash
vcftools --gzvcf input.vcf.gz --het
```

### Inbreeding

BCFtools v1.15.1. ROH classified as short (100 kb – 1 Mb) or long (> 1 Mb). F_ROH = sum of ROH lengths / total autosomal genome length (1,092,342 kb).

```bash
bcftools roh --rec-rate 2e-8 --AF-dflt 0.4 -O r input.vcf.gz > roh_output.txt
```

### Contemporary Ne

ONeSAMP v3.0.1. Run on 10 largest chromosomes (~816–819 SNPs per population). VCF converted to Genepop with vcf2popgen. Mutation rate μ = 8.11 × 10⁻⁸ (Wang et al. 2019).

```bash
vcftools --gzvcf pop.vcf.gz --thin 1000000 --max-missing 1 --recode
```

### CSS (Ireland + Britain vs Sweden)

SNP-by-SNP F_ST:

```bash
vcftools --gzvcf input.vcf.gz --weir-fst-pop IreBrit.txt --weir-fst-pop Swe.txt
```

Allele frequencies:

```bash
vcftools --gzvcf input.vcf.gz --freq --keep pop.txt
```

XP-EHH (selscan v1.2.0, per chromosome):

```bash
selscan --xpehh --vcf irebrit_chrN.vcf --vcf-ref swe_chrN.vcf --map chrN.map
norm --xpehh --files chrN.xpehh.out --bins 100 --bp-win --winsize 50000 --min-snps 10
```

### pcadapt

v4.4.1, qvalue v2.34.0. K = 2, FDR 0.01. Outlier positions extracted from `.bim` file, intersected with annotation via `extract_genes_pcadapt.py`. See `pcadapt_outliers.R`.

### GO overrepresentation

Pooled outlier gene lists submitted to g:Profiler. Query species: chicken (*Gallus gallus*). Significance threshold: g:SCS-adjusted P_adj < 0.05.
