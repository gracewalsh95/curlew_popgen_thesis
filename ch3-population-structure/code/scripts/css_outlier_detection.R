# css_outlier_detection.R
#
# Identifies CSS outlier SNPs, validates with flanking support,
# defines genomic regions, and merges nearby regions.
#
# Input:  CSS_smoothed_20kb.txt (output of smoothing_window.py)
# Output: selection_SNPs.txt, merged_regions.txt
#
# Key parameters:
#   - Core threshold: top 0.1% of smoothed CSS scores
#   - Support threshold: top 1% of smoothed CSS scores
#   - Flanking distance: 500 kb (each SNP must have >= 5 supporting SNPs within this distance)
#   - Region merge distance: 200 kb

library(data.table)
library(tidyverse)

# Parameters
CORE_QUANTILE   <- 0.999    # Top 0.1%
SUPPORT_QUANTILE <- 0.99    # Top 1%
FLANK_DIST      <- 500000   # 500 kb
MIN_SUPPORT     <- 5        # Minimum flanking SNPs in top 1%
MERGE_DIST      <- 200000   # Merge regions within 200 kb

# Read smoothed CSS data
data <- fread("CSS_smoothed_20kb.txt", header = TRUE)
data <- data[order(data$chromosome, data$position), ]

# Define thresholds
thresh_core <- quantile(data$CSS_smoothed, CORE_QUANTILE, na.rm = TRUE)
thresh_supp <- quantile(data$CSS_smoothed, SUPPORT_QUANTILE, na.rm = TRUE)

candidates <- data %>% filter(CSS_smoothed > thresh_core)
supporting <- data %>% filter(CSS_smoothed > thresh_supp)

# Validate: each core SNP must have >= MIN_SUPPORT supporting SNPs within FLANK_DIST ---
validated_snps <- candidates %>%
  mutate(has_support = map_lgl(seq_len(n()), function(i) {
    support_count <- supporting %>%
      filter(
        chromosome == chromosome[i],
        position >= (position[i] - FLANK_DIST),
        position <= (position[i] + FLANK_DIST),
        position != position[i]
      ) %>%
      nrow()
    support_count >= MIN_SUPPORT
  })) %>%
  filter(has_support == TRUE) %>%
  mutate(id = paste0(chromosome, ":", position)) %>%
  select(chromosome, position, id, CSS_smoothed)

cat("Validated SNPs:", nrow(validated_snps), "\n")
fwrite(validated_snps, "selection_SNPs.txt", sep = "\t")

# Define regions from outermost flanking SNPs 
all_regions <- data.frame()

if (nrow(validated_snps) > 0) {
  for (i in seq_len(nrow(validated_snps))) {
    chr <- validated_snps$chromosome[i]
    pos <- validated_snps$position[i]

    flanking_snps <- supporting %>%
      filter(chromosome == chr,
             position >= (pos - FLANK_DIST),
             position <= (pos + FLANK_DIST))

    all_regions <- rbind(all_regions, data.frame(
      chromosome   = chr,
      start        = min(flanking_snps$position),
      end          = max(flanking_snps$position),
      core_snp_pos = pos,
      n_flanking   = nrow(flanking_snps)
    ))
  }
}

# Merge nearby regions
all_regions <- all_regions %>% arrange(chromosome, start)
merged_regions <- data.frame()

for (chr in unique(all_regions$chromosome)) {
  chr_regions <- all_regions %>% filter(chromosome == chr)
  current_start <- chr_regions$start[1]
  current_end   <- chr_regions$end[1]

  if (nrow(chr_regions) > 1) {
    for (i in 2:nrow(chr_regions)) {
      if (chr_regions$start[i] - current_end <= MERGE_DIST) {
        current_end <- max(current_end, chr_regions$end[i])
      } else {
        merged_regions <- rbind(merged_regions, data.frame(
          chromosome = chr, start = current_start, end = current_end,
          region_size = current_end - current_start
        ))
        current_start <- chr_regions$start[i]
        current_end   <- chr_regions$end[i]
      }
    }
  }

  merged_regions <- rbind(merged_regions, data.frame(
    chromosome = chr, start = current_start, end = current_end,
    region_size = current_end - current_start
  ))
}

fwrite(merged_regions, "merged_regions.txt", sep = "\t")

cat("Initial regions:", nrow(all_regions), "\n")
cat("Merged regions:", nrow(merged_regions), "\n")
cat("Median region size:", median(merged_regions$region_size), "\n")
