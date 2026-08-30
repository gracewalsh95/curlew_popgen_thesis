# tajima_d_jackknife.R
#
# Tests whether mean Tajima's D differs significantly between populations
# using a block jackknife approach with non-overlapping 1 Mb blocks.
#
# Input:  combined_pixy_tajima_d.txt (concatenated pixy Tajima's D output)
# Output: Per-population estimates (mean, SE, 95% CI) and pairwise p-values
#
# Method:
#   1. Assign pixy windows to 1 Mb non-overlapping blocks (per chromosome)
#   2. Calculate mean Tajima's D per block
#   3. Estimate genome-wide mean and jackknife SE from block means
#   4. Pairwise differences tested with a two-sample t-statistic using
#      jackknife SEs combined(no equal-variance assumption)

library(tidyverse)

# Read pixy output
indtd <- read.table("combined_pixy_tajima_d.txt", sep = "\t", header = TRUE)

# Filter
indtd_filtered <- indtd %>%
  filter(!is.na(tajima_d), no_sites > 0)

# --- Assign 1 Mb blocks ---
indtd_filtered <- indtd_filtered %>%
  mutate(block = paste(chromosome, floor(window_pos_1 / 1e6), sep = "_"))

# Block jackknife function (operates on block means)
block_jackknife <- function(df) {
  block_means <- df %>%
    group_by(block) %>%
    summarise(block_mean = mean(tajima_d, na.rm = TRUE), .groups = "drop")

  n <- nrow(block_means)
  overall_mean <- mean(block_means$block_mean)

  # Leave-one-block-out
  jack_means <- sapply(1:n, function(i) {
    mean(block_means$block_mean[-i])
  })

  jack_se <- sqrt((n - 1) / n * sum((jack_means - mean(jack_means))^2))
  jack_ci <- 1.96 * jack_se

  return(data.frame(mean = overall_mean, se = jack_se, ci95 = jack_ci, n_blocks = n))
}

# Apply per population
results <- indtd_filtered %>%
  group_by(pop) %>%
  group_modify(~ block_jackknife(.x))

print(results)

# Pairwise t-tests
pops <- unique(indtd_filtered$pop)
cat("\nPairwise comparisons (block-mean jackknife):\n")
for (i in 1:(length(pops) - 1)) {
  for (j in (i + 1):length(pops)) {
    pop1 <- results[results$pop == pops[i], ]
    pop2 <- results[results$pop == pops[j], ]
    t_stat <- (pop1$mean - pop2$mean) / sqrt(pop1$se^2 + pop2$se^2)
    df <- min(pop1$n_blocks, pop2$n_blocks) - 1
    p_val <- 2 * pt(-abs(t_stat), df = df)
    cat(pops[i], "vs", pops[j], ": t =", round(t_stat, 3),
        ", df =", df, ", p =", format(p_val, digits = 4), "\n")
  }
}
