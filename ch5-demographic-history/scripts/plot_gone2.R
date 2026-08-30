#!/usr/bin/env Rscript
# plot_gone2.R — GONE2 Ne trajectories from subsampled replicates
# Two-panel figure: Britain + Ireland (panmictic) and Sweden (-x).
# Harmonic mean of Ne across 10 replicates per generation, with
# 2.5th–97.5th percentile envelope (95% CI).

library(ggplot2)
library(dplyr)
library(patchwork)

source("plot_theme.R")
# plot_theme.R defines: theme_ch4 (and palette used by other figures)

# paths (set to your working directories)
dir_bi  <- "path/to/britain_ireland_subsample_reps"
dir_swe <- "path/to/sweden_subsample_reps"
out     <- "path/to/output"

# constants
y_lab      <- expression(italic(N)[e] ~ "(×10"^3 * ")")
gen_time   <- 9.5
present    <- 2023
anchor_yr  <- 1500
anchor_gen <- (present - anchor_yr) / gen_time   # ~55 generations

# reader for -x mixed-format output
# GONE2 -x output has a commented header and 5 columns:
#   Rec_rate_bin, generation, N_T, Ne_metapop, d2
# Only generation (col 2) and Ne_metapop (col 4) are used.
read_mix <- function(path) {
  d <- read.table(path, header = FALSE, fill = TRUE, comment.char = "#")
  d <- d[!is.na(suppressWarnings(as.numeric(d[[5]]))), ]
  data.frame(Generation = as.numeric(d[[2]]), Ne = as.numeric(d[[4]]))
}

## Britain + Ireland
files_bi <- list.files(dir_bi, pattern = "^sub_rep.*_GONE2_Ne$", full.names = TRUE)
reps_bi  <- bind_rows(lapply(files_bi, function(f) {
  d <- read.table(f, header = TRUE); names(d)[1:2] <- c("Generation", "Ne"); d$rep <- f; d
}))
d_bi <- reps_bi %>%
  group_by(Generation) %>%
  summarise(Ne_hm = length(Ne) / sum(1 / Ne),
            lo = quantile(Ne, 0.025), hi = quantile(Ne, 0.975), .groups = "drop") %>%
  filter(Generation >= 5, Generation <= 150)

## Sweden
files_swe <- list.files(dir_swe, pattern = "^sub_rep.*_GONE2_Ne_mix$", full.names = TRUE)
reps_swe  <- bind_rows(lapply(files_swe, function(f) {
  d <- read_mix(f)
  d$rep <- f
  d
}))
d_swe <- reps_swe %>%
  group_by(Generation) %>%
  summarise(Ne_hm = length(Ne) / sum(1 / Ne),
            lo = quantile(Ne, 0.025), hi = quantile(Ne, 0.975), .groups = "drop") %>%
  filter(Generation >= 5, Generation <= 150)

## panels
p_bi <- ggplot(d_bi, aes(Generation, Ne_hm)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = anchor_gen, linetype = "dashed", colour = "grey50") +
  annotate("text", x = anchor_gen, y = Inf, label = anchor_yr,
           vjust = 1.4, hjust = -0.15, size = 3, colour = "grey50") +
  theme_ch4 +
  labs(x = "Generations ago", y = y_lab) +
  scale_y_continuous(labels = function(x) x / 1000,
                     expand = expansion(mult = c(0.02, 0.05))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  theme(axis.title.y = element_text(size = rel(0.9), margin = margin(r = 2)),
        axis.title.x = element_text(margin = margin(t = 2)),
        plot.margin  = margin(4, 6, 4, 4))

p_swe <- ggplot(d_swe, aes(Generation, Ne_hm)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = anchor_gen, linetype = "dashed", colour = "grey50") +
  annotate("text", x = anchor_gen, y = Inf, label = anchor_yr,
           vjust = 1.4, hjust = -0.15, size = 3, colour = "grey50") +
  theme_ch4 +
  labs(x = "Generations ago", y = y_lab) +
  scale_y_continuous(labels = function(x) x / 1000,
                     expand = expansion(mult = c(0.02, 0.05))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  theme(axis.title.y = element_text(size = rel(0.9), margin = margin(r = 2)),
        axis.title.x = element_text(margin = margin(t = 2)),
        plot.margin  = margin(4, 6, 4, 4))

combined <- p_bi + p_swe + plot_annotation(tag_levels = "a")
combined

## export
ggsave(file.path(out, "GONE2_IB_Sweden.png"), combined, width = 12, height = 6, dpi = 600)
ggsave(file.path(out, "GONE2_IB_Sweden.svg"), combined, width = 8, height = 4)
