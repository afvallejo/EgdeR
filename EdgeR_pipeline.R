#!/usr/bin/env Rscript

suppressMessages({
  library(edgeR)
  library(RColorBrewer)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript EdgeR_pipeline.R <counts_csv> <metadata_csv> <sample_column> <group_column>")
}
counts_path <- args[1]
metadata_path <- args[2]
sample_col <- args[3]
group_col <- args[4]

counts <- read.csv(counts_path, row.names = 1, check.names = FALSE)
metadata <- read.csv(metadata_path, stringsAsFactors = FALSE)

samples <- metadata[[sample_col]]
group <- factor(metadata[[group_col]])

if (length(group) != ncol(counts)) {
  stop("Metadata and counts dimensions do not match")
}

colnames(counts) <- samples

dge <- DGEList(counts = counts, group = group)
keep <- filterByExpr(dge, group = group)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

design <- model.matrix(~0 + group)
dge <- estimateDisp(dge, design)
fit <- glmFit(dge, design)

contrast <- if (ncol(design) == 2) c(1, -1) else rep(1, ncol(design))
lrt <- glmLRT(fit, contrast = contrast)
top <- topTags(lrt)
write.csv(top$table, file = "differential_expression.csv")

log_cpm <- cpm(dge, log = TRUE)
cols <- brewer.pal(max(3, nlevels(group)), "Set1")[group]

pdf("boxplot_logCPM.pdf")
boxplot(log_cpm, las = 2, col = cols, main = "logCPM")
dev.off()

pca <- prcomp(t(log_cpm))
pdf("PCA.pdf")
plot(pca$x[,1], pca$x[,2], pch = 19, col = cols,
     xlab = "PC1", ylab = "PC2", main = "PCA")
legend("bottomright", legend = levels(group), col = brewer.pal(max(3, nlevels(group)), "Set1"), pch = 19)
dev.off()

print(head(top$table))
