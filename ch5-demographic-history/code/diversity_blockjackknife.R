#!/usr/bin/env Rscript
# diversity_blockjackknife.R — genome-wide pi, Watterson's theta, Tajima's D
# Point estimates with 95% CIs from 1 Mb delete-one block jackknife.
# Pi: ratio estimator (sum diffs / sum comparisons)
# Theta_W: weighted mean by no_sites
# Tajima's D: mean of block means (unweighted — D is already normalised)
# Ne derived from theta_W / (4 * mu).

library(ggplot2)
library(dplyr)
library(patchwork)
library(ggsignif)

source("plot_theme.R")
# uses: theme_ch4, pop_colors, pop_order

# paths (set to your working directories)
data_dir <- "path/to/pixy_output"
out      <- "path/to/output"

pt_theme <- theme_ch4 + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# read
pi_raw    <- read.table(file.path(data_dir, "combined_pixy_pi.txt"),              sep = "\t", header = TRUE)
theta_raw <- read.table(file.path(data_dir, "combined_pixy_watterson_theta.txt"), sep = "\t", header = TRUE)
tajd_raw  <- read.table(file.path(data_dir, "combined_pixy_tajima_d.txt"),        sep = "\t", header = TRUE)

# block jackknife functions (1 Mb blocks)

# Pi: ratio estimator
bjk_pi <- function(df) {
  df <- df %>% 
    filter(!is.na(count_diffs), !is.na(count_comparisons), count_comparisons > 0) %>%
    mutate(block = paste(chromosome, floor(window_pos_1 / 1e6), sep = "_"))

  blocks <- unique(df$block)
  n <- length(blocks)
  overall <- sum(df$count_diffs) / sum(df$count_comparisons)

  jack_vals <- sapply(blocks, function(blk) {
    sub <- df[df$block != blk, ]
    sum(sub$count_diffs) / sum(sub$count_comparisons)
  })
  jack_se <- sqrt((n - 1) / n * sum((jack_vals - mean(jack_vals))^2))

  data.frame(mean_val = overall, se = jack_se,
             ci_lo = overall - 1.96 * jack_se,
             ci_hi = overall + 1.96 * jack_se,
             n_blocks = n)
}

# Theta_W: weighted mean by no_sites
bjk_theta <- function(df) {
  df <- df %>%
    filter(!is.na(avg_watterson_theta), no_sites > 0) %>%
    mutate(block = paste(chromosome, floor(window_pos_1 / 1e6), sep = "_"))

  blocks <- unique(df$block)
  n <- length(blocks)
  overall <- weighted.mean(df$avg_watterson_theta, df$no_sites)

  jack_vals <- sapply(blocks, function(blk) {
    sub <- df[df$block != blk, ]
    weighted.mean(sub$avg_watterson_theta, sub$no_sites)
  })
  jack_se <- sqrt((n - 1) / n * sum((jack_vals - mean(jack_vals))^2))

  data.frame(mean_val = overall, se = jack_se,
             ci_lo = overall - 1.96 * jack_se,
             ci_hi = overall + 1.96 * jack_se,
             n_blocks = n)
}

# Tajima's D: mean of block means
bjk_tajd <- function(df) {
  df <- df %>%
    filter(!is.na(tajima_d), no_sites > 0) %>%
    mutate(block = paste(chromosome, floor(window_pos_1 / 1e6), sep = "_"))

  block_means <- df %>%
    group_by(block) %>%
    summarise(bm = mean(tajima_d, na.rm = TRUE), .groups = "drop")

  n <- nrow(block_means)
  overall <- mean(block_means$bm)
  jack_vals <- sapply(1:n, function(i) mean(block_means$bm[-i]))
  jack_se <- sqrt((n - 1) / n * sum((jack_vals - mean(jack_vals))^2))

  data.frame(mean_val = overall, se = jack_se,
             ci_lo = overall - 1.96 * jack_se,
             ci_hi = overall + 1.96 * jack_se,
             n_blocks = n)
}

# compute estimates
pi_s    <- pi_raw    %>% group_by(pop) %>% group_modify(~ bjk_pi(.x))    %>% ungroup() %>% mutate(pop = factor(pop, levels = pop_order))
theta_s <- theta_raw %>% group_by(pop) %>% group_modify(~ bjk_theta(.x)) %>% ungroup() %>% mutate(pop = factor(pop, levels = pop_order))
tajd_s  <- tajd_raw  %>% group_by(pop) %>% group_modify(~ bjk_tajd(.x))  %>% ungroup() %>% mutate(pop = factor(pop, levels = pop_order))

# pairwise significance tests (two-sample t-test on jackknife SEs)
pairwise_bjk_test <- function(results_df) {
  pops <- levels(results_df$pop)
  out_list <- list()
  for (i in 1:(length(pops) - 1)) {
    for (j in (i + 1):length(pops)) {
      p1 <- results_df[results_df$pop == pops[i], ]
      p2 <- results_df[results_df$pop == pops[j], ]
      t_stat <- (p1$mean_val - p2$mean_val) / sqrt(p1$se^2 + p2$se^2)
      df <- min(p1$n_blocks, p2$n_blocks) - 1
      p_val <- 2 * pt(-abs(t_stat), df = df)
      out_list[[paste(pops[i], pops[j], sep = "_")]] <- data.frame(
        pop1 = pops[i], pop2 = pops[j],
        t = round(t_stat, 3), p = p_val,
        stringsAsFactors = FALSE
      )
    }
  }
  bind_rows(out_list)
}

pi_pw    <- pairwise_bjk_test(pi_s)
theta_pw <- pairwise_bjk_test(theta_s)
tajd_pw  <- pairwise_bjk_test(tajd_s)

cat("\n── Pi pairwise ──\n");      print(pi_pw)
cat("\n── Theta_W pairwise ──\n"); print(theta_pw)
cat("\n── Tajima's D pairwise ──\n"); print(tajd_pw)

# helper: p-value to annotation label
p_to_label <- function(p) {
  formatC(p, format = "e", digits = 1)
}

# panels with significance brackets
pt_panel <- function(df, ylab, pw, zeroline = FALSE,
                     show_pairs = list(c("Ireland", "Museum"), c("Britain", "Museum"))) {
  p <- ggplot(df, aes(pop, mean_val, colour = pop)) +
    { if (zeroline) geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") } +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2, linewidth = 0.7) +
    geom_point(size = 3) +
    scale_colour_manual(values = pop_colors) +
    labs(x = NULL, y = ylab) +
    pt_theme + theme(legend.position = "none")

  y_range <- max(df$ci_hi) - min(df$ci_lo)
  y_start <- max(df$ci_hi) + y_range * 0.05

  for (k in seq_along(show_pairs)) {
    pair <- show_pairs[[k]]
    row <- pw %>% filter((pop1 == pair[1] & pop2 == pair[2]) |
                         (pop1 == pair[2] & pop2 == pair[1]))
    if (nrow(row) == 1) {
      p <- p + geom_signif(
        comparisons = list(pair),
        annotations = p_to_label(row$p),
        y_position  = y_start + y_range * 0.08 * (k - 1),
        tip_length  = 0.01,
        textsize    = 3,
        colour      = "black"
      )
    }
  }
  p
}

p_pi    <- pt_panel(pi_s,    expression(pi),                    pi_pw)
p_theta <- pt_panel(theta_s, expression(theta[" W"]),           theta_pw)
p_d     <- pt_panel(tajd_s,  expression("Tajima's"~italic(D)),  pw = NULL, zeroline = TRUE, show_pairs = list())

trio <- p_pi + p_theta + p_d +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "plain"))
trio

ggsave(file.path(out, "diversity_pointCI_blockjacknife.png"), trio, width = 10, height = 4, dpi = 600)
ggsave(file.path(out, "diversity_pointCI_blockjacknife.pdf"), trio, width = 10, height = 4)

# summary table
summary_tbl <- bind_rows(
  pi_s    %>% mutate(metric = "pi"),
  theta_s %>% mutate(metric = "theta_w"),
  tajd_s  %>% mutate(metric = "tajima_d")
) %>%
  mutate(ci95 = (ci_hi - ci_lo) / 2) %>%
  select(metric, pop, mean_val, se, ci_lo, ci_hi, ci95, n_blocks) %>%
  arrange(metric, pop)

# Ne from theta_W
mu <- 8.11e-8
ne_tbl <- theta_s %>%
  transmute(metric = "Ne", pop,
            mean_val = mean_val / (4 * mu),
            se = se / (4 * mu),
            ci_lo = ci_lo / (4 * mu),
            ci_hi = ci_hi / (4 * mu),
            ci95 = (ci_hi - ci_lo) / 2,
            n_blocks = n_blocks)

summary_tbl <- bind_rows(summary_tbl, ne_tbl)

print(summary_tbl, n = Inf)
write.csv(summary_tbl, file.path(out, "diversity_summary_bjk.csv"), row.names = FALSE)

# pairwise test tables
pw_all <- bind_rows(
  pi_pw    %>% mutate(metric = "pi"),
  theta_pw %>% mutate(metric = "theta_w"),
  tajd_pw  %>% mutate(metric = "tajima_d")
)
write.csv(pw_all, file.path(out, "diversity_pairwise_bjk.csv"), row.names = FALSE)
