####################################################################
##  RNA-seq Normalisation Analysis with DESeq2 and edgeR
##  Using RNA_star_raw_count.rds data
####################################################################

### 0)  Install/load required packages -----------------------------------
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#     install.packages("BiocManager")
# }
# # Install required packages
# BiocManager::install(c("edgeR", "DESeq2"), ask = FALSE)
setwd("D:/bioinfor_normalization")
library(edgeR)
library(DESeq2)

### 1)  Load RNA_star data -----------------------------------------------
message("Loading RNA_star_raw_count.rds data...")
rds_file <- "RNA_star_raw_count.rds"

if (!file.exists(rds_file)) {
    stop("RNA_star_raw_count.rds file not found in current directory.")
}

cts <- readRDS(rds_file)

# Ensure it's a matrix
if (!is.matrix(cts)) {
    cts <- as.matrix(cts)
}

# Check if we need to clean up row names (remove any duplicates)
if (any(duplicated(rownames(cts)))) {
    rownames(cts) <- make.unique(rownames(cts))
}

### 2)  Data summary and quality check ------------------------------------
message("Count matrix loaded successfully!")
cat("Matrix dimension:", dim(cts)[1], "genes ×", dim(cts)[2], "samples\n\n")

### 3)  Determine filtering threshold based on data size -----------------
n_samples <- ncol(cts)
# Adaptive filtering: require expression in at least 10% of samples or minimum 2 samples
min_samples <- max(2, ceiling(n_samples * 0.1))

message("Applying filtering: CPM > 0.5 in at least ", min_samples, " samples")

### 4)  edgeR: TMM normalisation -----------------------------------------
message("Running edgeR TMM normalisation...")
y <- DGEList(cts)
print(dim(y))
# Filter genes with low counts
# keep <- rowSums(cpm(y) > 0.5) >= min_samples
# y <- y[keep, , keep.lib.sizes = FALSE]

y <- calcNormFactors(y) # checked that default is TMM normalization

edge_tab <- y$samples
edge_tab$mean_raw_per_gene <- edge_tab$lib.size / nrow(y)
edge_tab$invFactor <- 1 / edge_tab$norm.factors
edge_r <- cor(edge_tab$mean_raw_per_gene, edge_tab$invFactor)

### 5)  DESeq2: Median-of-ratios normalisation ---------------------------
message("Running DESeq2 median-of-ratios normalisation...")
# A design of ~1 is sufficient as we are only normalizing
dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData   = data.frame(row.names = colnames(cts)),
    design    = ~1
)
# Use the same filter as edgeR for a fair comparison
#dds <- dds[keep, ]
dds <- estimateSizeFactors(dds)
print(dim(dds))

deseq_tab <- data.frame(
    lib.size   = colSums(counts(dds)),
    sizeFactor = sizeFactors(dds)
)
deseq_tab$mean_raw_per_gene <- deseq_tab$lib.size / nrow(dds)
deseq_r <- cor(deseq_tab$mean_raw_per_gene, deseq_tab$sizeFactor)

### 6)  Summaries --------------------------------------------------------
cat("\n=== ANALYSIS SUMMARY ===\n")
cat("Data source: RNA_star_raw_count.rds\n")
cat("Original genes:", nrow(cts), "\n")
cat("Genes retained after filtering:", nrow(y), "\n")
cat("Samples:", ncol(cts), "\n\n")

cat("=== NORMALISATION RESULTS ===\n")
cat("Pearson correlation (mean counts per gene vs scaling factor):\n")
cat("  edgeR (mean vs 1/TMM) :", round(edge_r, 3), "\n")
cat("  DESeq2 (mean vs SF)   :", round(deseq_r, 3), "\n\n")

# Additional summary statistics
cat("=== SCALING FACTOR SUMMARIES ===\n")
cat("edgeR TMM factors:\n")
cat("  Range: [", round(min(y$samples$norm.factors), 3), ", ",
    round(max(y$samples$norm.factors), 3), "]\n",
    sep = ""
)
cat("  Mean: ", round(mean(y$samples$norm.factors), 3), "\n")
cat("  CV: ", round(sd(y$samples$norm.factors) / mean(y$samples$norm.factors), 3), "\n\n")

cat("DESeq2 size factors:\n")
cat("  Range: [", round(min(sizeFactors(dds)), 3), ", ",
    round(max(sizeFactors(dds)), 3), "]\n",
    sep = ""
)
cat("  Mean: ", round(mean(sizeFactors(dds)), 3), "\n")
cat("  CV: ", round(sd(sizeFactors(dds)) / mean(sizeFactors(dds)), 3), "\n\n")

### 7)  Scatter-plots ----------------------------------------------------
message("Generating comparison plots...")

op <- par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))

# edgeR plot
plot(edge_tab$mean_raw_per_gene, edge_tab$invFactor,
    pch = 19, col = "#3182bd90", cex = 0.8,
    xlab = "Mean Raw Counts per Gene",
    ylab = "1 / TMM Factor",
    main = sprintf("edgeR TMM Normalisation\nRNA_star (r = %.3f)", edge_r)
)
abline(h = 1, col = "grey70", lty = 2, lwd = 2)
grid(col = "lightgrey", lty = 3)

# DESeq2 plot
plot(deseq_tab$mean_raw_per_gene, deseq_tab$sizeFactor,
    pch = 19, col = "#de2d2690", cex = 0.8,
    xlab = "Mean Raw Counts per Gene",
    ylab = "DESeq2 Size Factor",
    main = sprintf("DESeq2 Median-of-Ratios\nRNA_star (r = %.3f)", deseq_r)
)
abline(h = 1, col = "grey70", lty = 2, lwd = 2)
grid(col = "lightgrey", lty = 3)

par(op)

### 8)  Optional: Save normalized counts ----------------------------------
cat("=== NORMALIZED COUNT EXPORT ===\n")
cat("To export normalized counts, uncomment the desired lines below:\n\n")

# Uncomment these lines to save normalized counts:
# edgeR_normalized <- cpm(y, normalized.lib.sizes = TRUE)
# deseq2_normalized <- counts(dds, normalized = TRUE)
# write.csv(edgeR_normalized, "edgeR_normalized_RNA_star.csv")
# write.csv(deseq2_normalized, "DESeq2_normalized_RNA_star.csv")

### 9)  Method comparison summary -----------------------------------------
cat("=== METHOD COMPARISON ===\n")
cat("Both methods aim to correct for differences in library sizes and composition.\n")
cat("Higher correlation indicates better relationship between library size and scaling factor.\n\n")

if (abs(edge_r) > abs(deseq_r)) {
    cat("edgeR TMM shows stronger correlation (|", round(edge_r, 3), "| > |", round(deseq_r, 3), "|)\n")
} else {
    cat("DESeq2 shows stronger correlation (|", round(deseq_r, 3), "| > |", round(edge_r, 3), "|)\n")
}

cat("\nAnalysis complete for RNA_star dataset!\n")
