## -----------------------------------------------------------
##  DESeq2 vignette – replicate the normalisation step only
##  Data set: "airway" (GSE52778; 8 libraries, 57 852 genes)
## -----------------------------------------------------------

## 1.  Install/load packages  ----
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install(c("DESeq2", "airway"), ask = FALSE)

library(DESeq2) # core method
library(airway) # provides the count matrix + metadata

## 2.  Load the airway SummarizedExperiment  ----
data(airway) # object is called `airway`
se <- airway
# quick check
dim(se) # 57852 genes x 8 samples
colData(se)[, c("cell", "dex")] # two factors used in the vignette

## 3.  Build a DESeqDataSet  ----
dds <- DESeqDataSet(se, design = ~dex)

## 4.  **Normalisation** : estimate size factors  ----
dds <- estimateSizeFactors(dds) # median-of-ratios (Anders & Huber 2010)

sizeFactors(dds) # show the 8 scaling factors

## 5.  Obtain normalised counts  ----
normCounts <- counts(dds, normalized = TRUE)

## Optional: write them out
# write.csv(normCounts, file = "airway_normalised_counts.csv")

## 6.  Sanity check – library size before/after  ----
rawTot <- colSums(counts(dds))
normTot <- colSums(normCounts)


res <- data.frame(
    rawTot      = rawTot,
    normTot     = normTot,
    sizeFactor  = sizeFactors(dds)
)

## derive mean library size using the TRUE number of genes
nGenes <- nrow(dds) # <- no more hard-coding
res$mean_lib.size <- res$rawTot / nGenes

## quick check / plot
print(res)
(rho <- cor(res$mean_lib.size, res$sizeFactor, method = "pearson"))

plot(res$mean_lib.size, res$sizeFactor,
    xlab = "mean raw counts per gene",
    ylab = "size factor",
    pch  = 19
)
abline(h = 1, col = "grey70", lty = 2)
legend("topleft",
    bty = "n",
    legend = sprintf("Pearson r = %.3f", rho)
)
