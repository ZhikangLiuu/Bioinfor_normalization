############################################################
##  TMM normalisation only  –  GSE60450 example
############################################################

## 1.  Load edgeR ---------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("edgeR", ask = FALSE)

library(edgeR)

## 2.  Download the public gene-level count matrix ------------------------
fn <- "GSE60450_Lactation-GenewiseCounts.txt.gz"
if (!file.exists(fn)) {
    download.file(
        paste0(
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/",
            "GSE60450/suppl/", fn
        ),
        fn,
        mode = "wb"
    )
}

cts <- read.delim(fn, row.names = "EntrezGeneID")
cts <- cts[, -1] # drop gene-length column

## 3.  Build a DGEList ----------------------------------------------------
# sample labels in the same order as the publication
grp <- factor(c(
    "B.virgin", "B.virgin",
    "B.pregnant", "B.pregnant",
    "B.lactating", "B.lactating",
    "L.virgin", "L.virgin",
    "L.pregnant", "L.pregnant",
    "L.lactating", "L.lactating"
))

y <- DGEList(cts, group = grp)

## 4.  (Optional) light filtering – keeps ~60 % of genes ------------------
keep <- rowSums(cpm(y) > 0.5) >= 2
y <- y[keep, , keep.lib.sizes = FALSE]

## 5.  TMM normalisation --------------------------------------------------
y <- calcNormFactors(y) # the whole point of this script!

## 6.  Inspect the scaling factors ---------------------------------------
y$samples
y$samples$eff.libsize <- y$samples$lib.size * y$samples$norm.factors
y$samples


res <- y$samples # lib.size and norm.factors already here

nGenes <- nrow(y) # TRUE gene count → no hard code
res$mean_lib.size <- res$lib.size / nGenes
res$invNormFactor <- 1 / res$norm.factors # optional – often easier to eyeball

## 2) correlation ---------------------------------------------------------
(rho <- cor(res$mean_lib.size, res$invNormFactor, method = "pearson"))
# or cor.test(res$mean_lib.size, res$invNormFactor)

## 3) plot ----------------------------------------------------------------
plot(res$mean_lib.size, res$invNormFactor,
    xlab = "mean raw counts per gene",
    ylab = "1 / TMM norm.factor",
    pch = 19, col = "steelblue"
)
abline(h = 1, col = "grey70", lty = 2)
legend("topright",
    bty = "n",
    legend = sprintf("Pearson r = %.3f", rho)
)
