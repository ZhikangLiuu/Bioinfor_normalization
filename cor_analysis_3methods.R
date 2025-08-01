####################################################################
##  Normalisation comparison: edgeR-TMM, DESeq2, edgeR-UpperQuartile
##  Data source : RNA_star_raw_count.rds
####################################################################

### 0)  Load packages -------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("edgeR", "DESeq2"), ask = FALSE)

library(edgeR)
library(DESeq2)
library(limma)

### 1)  Read count matrix ---------------------------------------------------
cts_file <- "RNA_star_raw_count.rds"
message("Reading ", cts_file, " …")
if (!file.exists(cts_file)) stop("File not found: ", cts_file)
cts <- readRDS(cts_file)
if (!is.matrix(cts)) cts <- as.matrix(cts)

# unique gene identifiers ---------------------------------------------------
if (any(duplicated(rownames(cts)))) {
    rownames(cts) <- make.unique(rownames(cts))
}

cat("Matrix dimension:", nrow(cts), "genes ×", ncol(cts), "samples\n\n")

### 2)  edgeR object + filtering -------------------------------------------
y_raw <- DGEList(cts)

# filter: CPM > 0.5 in ≥ 10 % of samples (≥2)
min_samples <- max(2, ceiling(0.1 * ncol(cts)))
keep <- rowSums(cpm(y_raw) > 0.5) >= min_samples
y_raw <- y_raw[keep, , keep.lib.sizes = FALSE]

cat("Genes retained after filtering:", nrow(y_raw), "\n\n")

### 3)  edgeR TMM normalisation -------------------------------------------
y_tmm <- y_raw
y_tmm <- calcNormFactors(y_tmm, method = "TMM")

tmm_tab <- y_tmm$samples
tmm_tab$mean_raw_per_gene <- tmm_tab$lib.size / nrow(y_tmm)
tmm_tab$invFactor <- 1 / tmm_tab$norm.factors

### 4)  voom transformation (TMM + logCPM) -------------------------------
# use intercept-only design because we are only interested in normalisation
voom_design <- matrix(1, ncol(cts), 1)
colnames(voom_design) <- "Intercept"

v <- voom(y_tmm, voom_design, plot = FALSE) # uses the TMM factors from y_tmm

voom_tab <- data.frame(
    lib.size            = tmm_tab$lib.size,
    mean_raw_per_gene   = tmm_tab$mean_raw_per_gene,
    mean_voom_weight    = colMeans(v$weights)
)

### 5)  DESeq2 median-of-ratios -------------------------------------------
dds <- DESeqDataSetFromMatrix(
    countData = y_raw$counts,
    colData   = data.frame(row.names = colnames(y_raw$counts)),
    design    = ~1
)

dds <- estimateSizeFactors(dds)

deseq_tab <- data.frame(
    lib.size   = colSums(counts(dds)),
    sizeFactor = sizeFactors(dds)
)
deseq_tab$mean_raw_per_gene <- deseq_tab$lib.size / nrow(dds)

### 6)  Correlations --------------------------------------------------------
cor_tmm <- cor(tmm_tab$mean_raw_per_gene, tmm_tab$invFactor)
cor_voom <- cor(voom_tab$mean_raw_per_gene, voom_tab$mean_voom_weight)
cor_deseq <- cor(deseq_tab$mean_raw_per_gene, deseq_tab$sizeFactor)

cat("Pearson correlations (mean counts per gene vs scaling factor / weights)\n")
cat("  edgeR-TMM :", round(cor_tmm, 3), "\n")
cat("  voom     :", round(cor_voom, 3), "\n")
cat("  DESeq2   :", round(cor_deseq, 3), "\n\n")

### 7)  Plots ---------------------------------------------------------------
op <- par(mfrow = c(1, 3), mar = c(5, 5, 3, 2))

# edgeR TMM
plot(tmm_tab$mean_raw_per_gene, tmm_tab$invFactor,
    pch = 19, col = "#3182bd90", cex = 0.8,
    xlab = "Mean Raw Counts per Gene", ylab = "1 / TMM factor",
    main = sprintf("edgeR-TMM\nr = %.2f", cor_tmm)
)
abline(h = 1, col = "grey70", lty = 2)

# voom weights
plot(voom_tab$mean_raw_per_gene, voom_tab$mean_voom_weight,
    pch = 19, col = "#756bb190", cex = 0.8,
    xlab = "Mean Raw Counts per Gene", ylab = "Mean voom weight",
    main = sprintf("voom (TMM + logCPM)\nr = %.2f", cor_voom)
)
abline(h = mean(voom_tab$mean_voom_weight), col = "grey70", lty = 2)

# DESeq2
plot(deseq_tab$mean_raw_per_gene, deseq_tab$sizeFactor,
    pch = 19, col = "#de2d2690", cex = 0.8,
    xlab = "Mean Raw Counts per Gene", ylab = "DESeq2 size factor",
    main = sprintf("DESeq2\nr = %.2f", cor_deseq)
)
abline(h = 1, col = "grey70", lty = 2)

par(op)


