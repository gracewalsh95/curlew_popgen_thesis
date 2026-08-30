# plot_venn.R — Venn diagrams of shared/private predicted deleterious variants
#
# Reads bcftools isec sites.txt and plots overlap across three populations.
# bcftools isec was run in order: Brit, Ire, Swe
# The binary CODE column is positional: position 1 = Brit, 2 = Ire, 3 = Swe
#
# Input: bcftools isec output directories with sites.txt files
#   ./intersection/           combined high + moderate
#   ./intersection/high/      high-impact only
#   ./intersection/mod/       moderate-impact only

library(ggVennDiagram)
library(ggplot2)
library(patchwork)

setwd("/path/to/working/directory")

# combined high + moderate
sites <- read.table("intersection/sites.txt", header = FALSE,
                    colClasses = c("character", "integer", "character", "character", "character"))
colnames(sites) <- c("CHROM", "POS", "REF", "ALT", "CODE")

brit <- paste0(sites$CHROM[grepl("^1", sites$CODE)], ":", sites$POS[grepl("^1", sites$CODE)])
ire  <- paste0(sites$CHROM[grepl("^.1", sites$CODE)], ":", sites$POS[grepl("^.1", sites$CODE)])
swe  <- paste0(sites$CHROM[grepl("1$", sites$CODE)], ":", sites$POS[grepl("1$", sites$CODE)])

plot_combined <- ggVennDiagram(list(Britain = brit, Ireland = ire, Sweden = swe), label_size = 3) +
  scale_fill_gradient(low = "white", high = "lightblue")
ggsave("venn_combined.png", plot_combined, width = 5, height = 5)

# high impact only
sites <- read.table("intersection/high/sites.txt", header = FALSE,
                    colClasses = c("character", "integer", "character", "character", "character"))
colnames(sites) <- c("CHROM", "POS", "REF", "ALT", "CODE")

brit <- paste0(sites$CHROM[grepl("^1", sites$CODE)], ":", sites$POS[grepl("^1", sites$CODE)])
ire  <- paste0(sites$CHROM[grepl("^.1", sites$CODE)], ":", sites$POS[grepl("^.1", sites$CODE)])
swe  <- paste0(sites$CHROM[grepl("1$", sites$CODE)], ":", sites$POS[grepl("1$", sites$CODE)])

plot_high <- ggVennDiagram(list(Britain = brit, Ireland = ire, Sweden = swe), label_size = 3) +
  scale_fill_gradient(low = "white", high = "lightblue")

# moderate impact only
sites <- read.table("intersection/mod/sites.txt", header = FALSE,
                    colClasses = c("character", "integer", "character", "character", "character"))
colnames(sites) <- c("CHROM", "POS", "REF", "ALT", "CODE")

brit <- paste0(sites$CHROM[grepl("^1", sites$CODE)], ":", sites$POS[grepl("^1", sites$CODE)])
ire  <- paste0(sites$CHROM[grepl("^.1", sites$CODE)], ":", sites$POS[grepl("^.1", sites$CODE)])
swe  <- paste0(sites$CHROM[grepl("1$", sites$CODE)], ":", sites$POS[grepl("1$", sites$CODE)])

plot_mod <- ggVennDiagram(list(Britain = brit, Ireland = ire, Sweden = swe), label_size = 3) +
  scale_fill_gradient(low = "white", high = "lightblue")

# combined panel
plot_mod + plot_high + plot_layout(ncol = 2) + plot_annotation(tag_levels = "a")
ggsave("venn_mod_high.png", width = 10, height = 5)
