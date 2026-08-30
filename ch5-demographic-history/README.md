
# Chapter 5: Temporal Genomic Erosion

Code and parameters to reproduce the results in Chapter 5. Standard bioinformatics pipeline steps (read trimming, mapping, variant calling, quality filtering) follow the same methods described in Chapter 3. The following analyses also use identical methods to previous chapters and are documented there:

- PCA, ADMIXTURE → Chapter 3 supplementary
- Heterozygosity (`vcftools --het`) → Chapter 3 supplementary
- ROH calling and F_ROH → Chapter 4 supplementary (`roh_analysis.R`)
- Genetic load (SnpEff) → Chapter 4 supplementary (`build_load_table.sh`)
- Rxy block jackknife → Chapter 4 supplementary (`rxy_jackknife.R`)

Chapter 5 uses an additional 7 museum samples.

**Reference genome:** GCA_964106895.1 (*Numenius arquata*)

**Samples:** 63 total (56 contemporary, 7 museum). Mus1 excluded from temporal comparisons (clusters with Sweden on PCA). Curl1 excluded from per-individual estimates (reference-match artefact). Final temporal comparison: Ireland n=23, Britain n=18, Museum n=6, Sweden n=15.

## Scripts

| Script | Description |
|--------|-------------|
| `plot_theme.R` | Shared ggplot2 theme, colour palette, and helper functions. Source before running other R scripts |
| `plot_gone2.R` | Plots GONE2 Ne trajectories from 10 subsampled replicates. Computes harmonic mean of Ne per generation and 95% CI (2.5th–97.5th percentile envelope). Includes parser for GONE2 `-x` mixed-format output (Sweden) |
| `diversity_blockjackknife.R` | Calculates genome-wide π, θ_W, and Tajima's D with 95% CIs from 1 Mb delete-one block jackknife. Ratio estimator for π, weighted mean for θ_W, mean of block means for Tajima's D. Pairwise two-sample t-tests on jackknife SEs. Derives Ne from θ_W / (4μ) |

## Analysis Parameters

### Museum sample DNA damage assessment and correction

mapDamage2 v2.3.0 (Jónsson et al. 2013) on duplicate-marked BAMs. USER enzyme treatment applied during library preparation. Post-rescaling δ_d values ranged from 0.0056 to 0.0075.

```bash
mapDamage -i sample.duplMarked.bam \
  -r reference.fasta \
  --rescale --merge-libraries
```

Mus6 downsampled to 78% of reads prior to rescaling:

```bash
samtools view -h -s 0.78 -b Mus6.duplMarked.bam -o Mus6_downsampled.duplMarked.bam
samtools index Mus6_downsampled.duplMarked.bam
```

### PSMC

v0.6.5 (Li & Durbin 2011). Run on 6 samples: 2 Ireland/Britain, 2 Sweden (highest coverage contemporary), 2 museum (excluded from final figures due to insufficient coverage). Scaling: generation time = 9.5 years (Bird et al. 2020), mutation rate = 8.11 × 10⁻⁸ per site per generation (Wang et al. 2019; Charadriiformes).

**Step 1:** Filter BAMs to autosomes, remove unmapped/secondary/duplicate/QC-fail reads.

```bash
samtools view -b -F 3844 -q 25 -L autosomes.bed input.bam -o sample_filtered.bam
samtools index sample_filtered.bam
```

**Step 2:** Assess per-sample autosomal coverage with Qualimap2 v2.3.

```bash
qualimap bamqc -bam filtered.bam -gff autosomes.bed -outdir outdir
```

**Step 3:** Generate diploid consensus FASTQ (min depth 10, max depth 2× mean coverage).

```bash
bcftools mpileup -C 50 -f reference.fasta filtered.bam | \
  bcftools call -c | \
  vcfutils.pl vcf2fq -d 10 -D 2x_mean_coverage | \
  gzip > sample_diploid.fq.gz
```

**Step 4:** Check proportion of missing/low-quality bases.

```bash
TOTAL=$(zcat sample.fq.gz | awk 'NR%4==2' | wc -c)
MISSING=$(zcat sample.fq.gz | awk 'NR%4==2' | tr -cd 'nN' | wc -c)
PCT=$(echo "scale=2; $MISSING/$TOTAL*100" | bc)
```

**Step 5:** Convert to PSMC input format.

```bash
fq2psmcfa -q 20 sample_diploid.fq.gz > sample.psmcfa
```

**Step 6:** Run PSMC.

```bash
psmc -N30 -t5 -r5 -p "4+30*2+4+6+10" -o sample.psmc sample.psmcfa
```

Artefact check following Hilgers et al. (2025):

```bash
psmc -N30 -t5 -r5 -p "2+2+30*2+4+6+10" -o sample_artefact_check.psmc sample.psmcfa
```

Parameters follow Bruniche-Olsen et al. (2021), applied across 68 avian species. Robustness confirmed across 12 parameter combinations on highest-coverage sample (Brit16, 25.48×).

**Step 7:** Bootstrap (100 replicates per sample).

```bash
splitfa sample.psmcfa > sample_split.psmcfa

seq 100 | xargs -I{} -P 8 psmc -N30 -t5 -r5 -b \
  -p "4+30*2+4+6+10" \
  -o sample_round-{}.psmc sample_split.psmcfa
```

**Step 8:** Combine and plot.

```bash
cat sample.psmc sample_round-*.psmc > sample_combined.psmc
psmc_plot.pl -R -u 8.11e-8 -g 9.5 output_prefix sample_combined.psmc
```

### GONE2

v2.0 (Santiago et al. 2025). Contemporary samples only, restricted to macrochromosomes (chr 1–19). Related individuals removed (KING kinship coefficient > 0.0884). Britain + Ireland treated as panmictic; Sweden run with `-x` (structure model).

**Step 1:** Convert VCF to PLINK PED format.

```bash
plink --vcf noMAF_SNPs.vcf.gz \
  --keep population_keep.txt \
  --geno 0_or_0.1 \
  --chr 1-19 \
  --recode --allow-extra-chr --chr-set 40 \
  --set-missing-var-ids @:# \
  --out output_prefix
```

**Step 2:** Run GONE2.

```bash
# Britain + Ireland (panmictic, default model)
gone2 -g 0 -r 2 -u 0.05 -t 8 input.ped

# Sweden (structure-corrected model)
gone2 -g 0 -r 2 -u 0.05 -t 8 -x input.ped
```

**Step 3:** Subsampling for confidence intervals (10 replicates, 80,000 SNPs each).

```bash
# Britain + Ireland
for i in $(seq 1 10); do
  gone2 -g 0 -r 2 -u 0.05 -t 8 -s 80000 -S $i -o out_rep${i} input.ped
done

# Sweden (with -x)
for i in $(seq 1 10); do
  gone2 -g 0 -r 2 -u 0.05 -t 8 -s 80000 -x -S $i -o out_rep${i} input.ped
done
```

CI method: harmonic mean of Ne across 10 replicates per generation, with 2.5th–97.5th percentile envelope. Harmonic mean used because Ne is inversely related to drift. See `plot_gone2.R`.

Key parameter justification: r = 2 cM/Mb — restricted to macrochromosomes (chr 1–19), which have lower recombination than microchromosomes. r = 4 sensitivity shown in supplementary. Sensitivity runs tested combinations of: MAC (0 or 3), missing data (0% or 10%), recombination rate (2 or 4 cM/Mb), and structure model (default vs `-x`).

### Genetic diversity

pixy v2.0.0.beta8 (Korunes & Samuk 2021). Per chromosome, all-sites VCF with no MAF filter. Same pipeline as Chapter 3, with museum samples included.

```bash
pixy --stats pi watterson_theta tajima_d \
  --vcf allsites_chrN.vcf.gz \
  --populations popfile_4pop.tsv \
  --window_size 15000 \
  --n_cores 1
```

Concatenated across chromosomes:

```bash
cat chr1/pixy_pi.txt > combined_pixy_pi.txt
for c in $(seq 2 39); do tail -n +2 chr${c}/pixy_pi.txt >> combined_pixy_pi.txt; done
```

Block jackknife and pairwise tests: see `diversity_blockjackknife.R`.

### Significance testing

R v4.3.3. Per-individual metrics (F_ROH, genetic load ratios) compared between populations using two-sided Wilcoxon rank-sum tests (same as Chapter 4). Population-level diversity statistics (π, θ_W) compared using two-sample t-tests on block-jackknife standard errors (see `diversity_blockjackknife.R`). Tajima's D interpreted via per-cohort CIs; no pairwise testing applied.


## Notes

1. PSMC parameters (`-N30 -t5 -r5 -p "4+30*2+4+6+10"`) follow Bruniche-Olsen et al. (2021), applied across 68 avian species.
2. GONE2 `-x` output uses a 5-column mixed format requiring special parsing. See `read_mix()` in `plot_gone2.R`.
