# roh_analysis.R — ROH processing and statistical tests
#
# Reads raw BCFtools/RoH output, categorises ROH into short (100 kb - 1 Mb)
# and long (>1 Mb), calculates FROH, and runs statistical tests across
# three Eurasian curlew populations (Ireland, Britain, Sweden).
#
# Input:
#   BCFtools/RoH output file (--output-type r), with # header lines
#   Columns: RG, Sample, Chromosome, Start, End, Length_bp, Num_markers, Quality
#
# Tests:
#   Kruskal-Wallis + pairwise Wilcoxon (BH correction) for FROH, NROH,
#   and mean ROH length across populations
#   Spearman rank correlations of per-window (1 Mb) ROH proportions
#
# Output:
#   roh_statistical_tests.csv
#   Spearman correlation results printed to console
#
# Autosomal length: 1,092,342,054 bp (GCA_964106895.1)
# Curl1 excluded (reference genome individual, reference bias)

library(data.table)
library(dplyr)

setwd("/path/to/working/directory")

autosome_kb <- 1092342

assign_population <- function(data) {
  data$Pop <- ifelse(grepl("Brit", data$Sample), "Britain",
              ifelse(grepl("Curl", data$Sample), "Ireland",
              ifelse(grepl("Swe",  data$Sample), "Sweden", "Other")))
  return(data)
}

# Read and process raw ROH
ROHall <- read.table("bcftools_output.txt",
                     header = FALSE, comment.char = "#", stringsAsFactors = FALSE)
colnames(ROHall) <- c("RG", "Sample", "Chromosome", "Start", "End",
                       "Length_bp", "Num_markers", "Quality")

ROHall <- ROHall[ROHall$Sample != "Curl1", ]
ROHall <- assign_population(ROHall)
ROHall$KB <- ROHall$Length_bp / 1000
ROHall <- ROHall %>%
  mutate(ROH_category = case_when(
    KB > 100 & KB < 1000 ~ "short",
    KB >= 1000           ~ "long"
  )) %>%
  filter(!is.na(ROH_category))

# Summaries
short_roh <- ROHall[ROHall$ROH_category == "short", ]
long_roh  <- ROHall[ROHall$ROH_category == "long", ]

# all ROH
total_length <- aggregate(KB ~ Sample, data = ROHall, FUN = sum)
total_count <- aggregate(KB ~ Sample, data = ROHall, FUN = length)
names(total_count)[2] <- "count"
summary_data_all <- merge(total_length, total_count, by = "Sample")
summary_data_all <- assign_population(summary_data_all)

# short ROH
total_length_short <- aggregate(KB ~ Sample, data = short_roh, FUN = sum)
total_count_short <- aggregate(KB ~ Sample, data = short_roh, FUN = length)
names(total_count_short)[2] <- "count"
summary_data_short <- merge(total_length_short, total_count_short, by = "Sample")
summary_data_short <- assign_population(summary_data_short)
summary_data_short$MB <- summary_data_short$KB / 1000

# long ROH
total_length_long <- aggregate(KB ~ Sample, data = long_roh, FUN = sum)
total_count_long <- aggregate(KB ~ Sample, data = long_roh, FUN = length)
names(total_count_long)[2] <- "count"
summary_data_long <- merge(total_length_long, total_count_long, by = "Sample")
summary_data_long <- assign_population(summary_data_long)
summary_data_long$MB <- summary_data_long$KB / 1000

# FROH
summary_data_short$FROH <- summary_data_short$KB / autosome_kb
summary_data_long$FROH  <- summary_data_long$KB / autosome_kb
summary_data_all$FROH   <- summary_data_all$KB / autosome_kb

# mean ROH lengths per individual
meanshort <- tapply(short_roh$KB, short_roh$Sample, mean)
meanshort <- data.frame(Sample = names(meanshort), Kb = as.numeric(meanshort))
meanshort$MB <- meanshort$Kb / 1000
meanshort <- assign_population(meanshort)

meanslong <- tapply(long_roh$KB, long_roh$Sample, mean)
meanslong <- data.frame(Sample = names(meanslong), Kb = as.numeric(meanslong))
meanslong$MB <- meanslong$Kb / 1000
meanslong <- assign_population(meanslong)

# statistical tests 
# FROH
kw_short <- kruskal.test(FROH ~ Pop, data = summary_data_short)
kw_long  <- kruskal.test(FROH ~ Pop, data = summary_data_long)
kw_all   <- kruskal.test(FROH ~ Pop, data = summary_data_all)
pw_short <- pairwise.wilcox.test(summary_data_short$FROH, summary_data_short$Pop, p.adjust.method = "BH")
pw_long  <- pairwise.wilcox.test(summary_data_long$FROH, summary_data_long$Pop, p.adjust.method = "BH")
pw_all   <- pairwise.wilcox.test(summary_data_all$FROH, summary_data_all$Pop, p.adjust.method = "BH")

# NROH
kw_nroh_short <- kruskal.test(count ~ Pop, data = summary_data_short)
kw_nroh_long  <- kruskal.test(count ~ Pop, data = summary_data_long)
pw_nroh_short <- pairwise.wilcox.test(summary_data_short$count, summary_data_short$Pop, p.adjust.method = "BH")
pw_nroh_long  <- pairwise.wilcox.test(summary_data_long$count, summary_data_long$Pop, p.adjust.method = "BH")

# Mean ROH length
kw_mean_short <- kruskal.test(MB ~ Pop, data = meanshort)
kw_mean_long  <- kruskal.test(MB ~ Pop, data = meanslong)
pw_mean_short <- pairwise.wilcox.test(meanshort$MB, meanshort$Pop, p.adjust.method = "BH")
pw_mean_long  <- pairwise.wilcox.test(meanslong$MB, meanslong$Pop, p.adjust.method = "BH")

# results 
fmt_p <- function(p) ifelse(p < 0.001, "< 0.001", round(p, 3))

results <- data.frame(
  Metric = c("FROH_short", "FROH_long", "FROH_total",
             "NROH_short", "NROH_long",
             "Mean_length_short", "Mean_length_long"),
  KW_chi_sq = c(kw_short$statistic, kw_long$statistic, kw_all$statistic,
                kw_nroh_short$statistic, kw_nroh_long$statistic,
                kw_mean_short$statistic, kw_mean_long$statistic),
  KW_p = c(kw_short$p.value, kw_long$p.value, kw_all$p.value,
           kw_nroh_short$p.value, kw_nroh_long$p.value,
           kw_mean_short$p.value, kw_mean_long$p.value),
  Ire_vs_Brit = c(pw_short$p.value["Britain","Ireland"],
                  pw_long$p.value["Britain","Ireland"],
                  pw_all$p.value["Britain","Ireland"],
                  pw_nroh_short$p.value["Britain","Ireland"],
                  pw_nroh_long$p.value["Britain","Ireland"],
                  pw_mean_short$p.value["Britain","Ireland"],
                  pw_mean_long$p.value["Britain","Ireland"]),
  Ire_vs_Swe = c(pw_short$p.value["Sweden","Ireland"],
                 pw_long$p.value["Sweden","Ireland"],
                 pw_all$p.value["Sweden","Ireland"],
                 pw_nroh_short$p.value["Sweden","Ireland"],
                 pw_nroh_long$p.value["Sweden","Ireland"],
                 pw_mean_short$p.value["Sweden","Ireland"],
                 pw_mean_long$p.value["Sweden","Ireland"]),
  Brit_vs_Swe = c(pw_short$p.value["Sweden","Britain"],
                  pw_long$p.value["Sweden","Britain"],
                  pw_all$p.value["Sweden","Britain"],
                  pw_nroh_short$p.value["Sweden","Britain"],
                  pw_nroh_long$p.value["Sweden","Britain"],
                  pw_mean_short$p.value["Sweden","Britain"],
                  pw_mean_long$p.value["Sweden","Britain"])
)
results$KW_p        <- fmt_p(results$KW_p)
results$Ire_vs_Brit <- fmt_p(results$Ire_vs_Brit)
results$Ire_vs_Swe  <- fmt_p(results$Ire_vs_Swe)
results$Brit_vs_Swe <- fmt_p(results$Brit_vs_Swe)

write.csv(results, "roh_statistical_tests.csv", row.names = FALSE)
print(results)

# Spearman correlations of per-window (1 Mb) ROH proportions
runs <- as.data.table(ROHall)
runs$Chromosome <- as.character(runs$Chromosome)

pop_n <- c(Britain = 18, Ireland = 22, Sweden = 15)
pops  <- names(pop_n)

# 1 Mb bins per chromosome
chrom_lengths <- runs[, .(max_pos = max(End)), by = Chromosome]
bins <- chrom_lengths[, .(bin = seq(0, floor(max_pos / 1e6))), by = Chromosome]
bins[, mid := bin * 1e6 + 5e5]

# per-population ROH proportion per window
prop_list <- lapply(pops, function(pop) {
  pop_runs <- runs[Pop == pop]
  props <- sapply(seq_len(nrow(bins)), function(i) {
    chr <- bins$Chromosome[i]; pos <- bins$mid[i]
    sum(pop_runs$Chromosome == chr & pop_runs$Start <= pos & pop_runs$End >= pos) / pop_n[pop]
  })
  data.table(Chromosome = bins$Chromosome, bin = bins$bin, prop = props)
})
names(prop_list) <- pops

wide <- merge(prop_list[["Ireland"]], prop_list[["Britain"]],
              by = c("Chromosome", "bin"), suffixes = c("_Ire", "_Brit"))
wide <- merge(wide, prop_list[["Sweden"]], by = c("Chromosome", "bin"))
setnames(wide, "prop", "prop_Swe")

hw_nonzero <- wide[prop_Ire > 0 | prop_Brit > 0 | prop_Swe > 0]

for (pair in list(c("prop_Ire","prop_Brit","Ireland","Britain"),
                  c("prop_Ire","prop_Swe","Ireland","Sweden"),
                  c("prop_Brit","prop_Swe","Britain","Sweden"))) {
  ct <- cor.test(hw_nonzero[[pair[1]]], hw_nonzero[[pair[2]]], method = "spearman")
  cat(pair[3], "vs", pair[4], ": rho =", round(ct$estimate, 3),
      ", p =", format.pval(ct$p.value, digits = 3), "\n")
}
