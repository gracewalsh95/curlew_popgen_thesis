# pcadapt_outliers.R
#
# Detects outlier SNPs associated with population structure using pcadapt.
#
# Input:  PLINK .bed file (full SNP dataset, not LD-pruned)
# Output: Outlier SNPs (row numbers referencing the .bim file)
#
# Key parameters:
#   - K = 2 PCs (selected from scree plot)
#   - FDR threshold = 0.01 (Storey q-value method)
#
# Outlier positions are extracted from the .bim file using the row indices
# and intersected with the genome annotation using extract_genes_pcadapt.py.

library(pcadapt)
library(qvalue)

# Read input
filename <- read.pcadapt("input.bed", type = "bed")

# Inspect scree plot to choose K
x_scree <- pcadapt(input = filename, K = 10)
plot(x_scree, option = "screeplot")

# Run pcadapt with K = 2
x <- pcadapt(filename, K = 2)

# Apply FDR correction
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.01
outliers <- which(qval < alpha)

cat("Outlier SNPs at FDR <", alpha, ":", length(outliers), "\n")

# Export outlier indices
# These are row numbers in the .bim file
write.table(outliers, file = "pcadapt_outliers_FDR0.01.txt",
            sep = "\t", row.names = FALSE, col.names = FALSE)
