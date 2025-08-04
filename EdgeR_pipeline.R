%%writefile /content/EdgeR_pipeline.R
#!/usr/bin/env Rscript
.libPaths(c("/content/R_libraries_EdgeR/content/R_libraries/", .libPaths()))
# Load CRAN packages
# Load required libraries
library(devtools)
library(tidyverse)
library(dplyr)
library(edgeR)
library(limma)
library(tximport)
library(biomaRt)
library(genefilter)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(EnsDb.Mmusculus.v79)
library(rhdf5)
library(RColorBrewer)
library(vioplot)
library(gplots)
library(ggrepel)
library(rtracklayer)
library(FactoMineR)
library(factoextra)
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


design <- model.matrix(~0+group) # assigning matrix model to the binary configuration (control vs test)

#design <- model.matrix(~donor+group) ---> this would instead seek to lay out matrix model to investigate donor effect on the different samples

#print(design)
colnames(design) # observing the generated matrix model


dge <- DGEList(counts= counts,genes=row.names(counts)) # produces a list of DEG using the dge() function
dge_norm <- calcNormFactors(dge,lib.size=T,method="TMM") # Uses TMM normalisation for the dge list
cpm.count <- cpm(dge_norm) # cacluating cpm for the normalised DGE data
keep <- filterByExpr(dge_norm,group=group) # edgeR fucntion that using compiting to automatically filter genes by expression level
dge_norm_filt<-dge_norm[keep,keep.lib.sizes=FALSE]
#keep = rowSums(cpm.count > as.numeric(2) ) >= as.numeric(3)
#dge_norm_filt <- dge_norm[keep,,keep.lib.sizes=FALSE] # -----> edgeR filtered analysis

# Convert raw counts to logCPMs
logcounts1 <- cpm(dge_norm_filt,log=TRUE)

# Save the plot as a PDF
pdf(file = "boxplot_logcounts1.pdf")
boxplot(logcounts1,
        xlab="Sample Groups",
        ylab="Log2 counts per million",
        las=2,
        col=col.cell,
        notch=TRUE,
        main="Boxplots of logCPMs")
par(cex.axis=0.5, las=1)
points(jitter(1:length(logcounts1)), logcounts1, col="darkgray", pch=16)
abline(h=median(logcounts1), col="blue", lwd=2)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
dev.off() # Close the PDF device

# Display the same plot in the R environment
boxplot(logcounts1,
        xlab="Samples",
        ylab="Log2 counts per million",
        las=2,
        col=col.cell,
        notch=TRUE,
        main="Boxplots of logCPMs (TMM normalized)")
par(mar=c(7, 4, 4, 2) + 0.1) # Adjusting margin for better visibility
points(jitter(1:length(logcounts1)), logcounts1, col="darkgray", pch=16)
abline(h=median(logcounts1), col="blue", lwd=2)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")


plot_name <- "SamplesIQRvsMedian"
file_name <- paste0(root, "_", plot_name, ".pdf")

pdf(file_name)
#  Make IQR_Median Plots
counts_filtered <- dge_norm_filt$counts
log_counts<-log(counts_filtered+1)
CPM<-cpm(dge_norm_filt)

IQR<-apply(CPM, 2, IQR)
Median<-apply(CPM, 2, median)
diff1<-mean(Median)-min(Median)
diff2<-max(Median)-mean(Median)
diff3<-mean(IQR)-min(IQR)
diff4<-max(IQR)-mean(IQR)

#  These plot settings work for this data and may need to be adjusted for your own data.
Xlim=c(mean(Median)-2*diff1,mean(Median)+2*diff2)
Ylim=c(mean(IQR)-2*diff3,mean(IQR)+2*diff4)


plot(Median, IQR, main="IQR vs. Median TMM normaliztion", type="n", xlim=Xlim,ylim=Ylim)
text(Median, IQR, labels=names(IQR))

#  Make boxes for StDev.
Median_mean<-mean(Median)
c_sd1_mean<-sd(Median)
c_sd2_mean<-2*sd(Median)
c_sd3_mean<-3*sd(Median)
IQR_mean<-mean(IQR)
c_sd1_IQR<-sd(IQR)
c_sd2_IQR<-2*sd(IQR)
c_sd3_IQR<-3*sd(IQR)

x0_c<-Median_mean-c_sd1_mean
y0_c<-IQR_mean-c_sd1_IQR
x1_c<-Median_mean+c_sd1_mean
y1_c<-IQR_mean+c_sd1_IQR

x0_c.2<-Median_mean-c_sd2_mean
y0_c.2<-IQR_mean-c_sd2_IQR
x1_c.2<-Median_mean+c_sd2_mean
y1_c.2<-IQR_mean+c_sd2_IQR

x0_c.3<-Median_mean-c_sd3_mean
y0_c.3<-IQR_mean-c_sd3_IQR
x1_c.3<-Median_mean+c_sd3_mean
y1_c.3<-IQR_mean+c_sd3_IQR

segments(x0_c,y0_c, x1=x1_c, y1=y0_c, col="blue")
segments(x0_c,y0_c, x1=x0_c, y1=y1_c, col="blue")
segments(x1_c,y0_c, x1=x1_c, y1=y1_c, col="blue")
segments(x0_c,y1_c, x1=x1_c, y1=y1_c, col="blue")

segments(x0_c.2,y0_c.2, x1=x1_c.2, y1=y0_c.2, col="red")
segments(x0_c.2,y0_c.2, x1=x0_c.2, y1=y1_c.2, col="red")
segments(x1_c.2,y0_c.2, x1=x1_c.2, y1=y1_c.2, col="red")
segments(x0_c.2,y1_c.2, x1=x1_c.2, y1=y1_c.2, col="red")

#  This portion is out of range so I removed it.
segments(x0_c.3,y0_c.3, x1=x1_c.3, y1=y0_c.3, col="green")
segments(x0_c.3,y0_c.3, x1=x0_c.3, y1=y1_c.3, col="green")
segments(x1_c.3,y0_c.3, x1=x1_c.3, y1=y1_c.3, col="green")
segments(x0_c.3,y1_c.3, x1=x1_c.3, y1=y1_c.3, col="green")

outliers2 <- which(IQR < IQR_mean - 2*sd(IQR) | IQR > IQR_mean + 2*sd(IQR))
points(Median[outliers2], IQR[outliers2], col="red", pch=19)

dev.off()

cpm <- cpm(dge_norm_filt)
#cpm<-cpm(v)

Y <- apply(cpm, 1, function(y) scale(y, center=TRUE, scale=FALSE))
s <- svd(Y)
#  Get Variance explained by the components
Var1<-s$d^2/sum(s$d^2)
Var2<-cumsum(s$d^2/sum(s$d^2))
Var<-cbind(Var1,Var2)
Var<-round(Var,3)
colnames(Var)<-c("Percent Variance","Cumulative % Variance")




#dge_norm<- estimateGLMTagwiseDisp(dge_norm, design)
dge_norm <- estimateDisp(dge_norm, design, robust=TRUE)
dge_norm.fit <- glmFit(dge_norm, design) # fitting the normalised DGE data, after calculating dispersion values
dge.lrt <- glmLRT(dge_norm.fit,contrast=c(1,-1)) # likelihood test, contrasting infected samples against healthy samples
#dge.lrt <- glmLRT(dge_norm.fit,contrast=c(1,-1,0,0)) # likelihood test, contrasting infected samples against healthy samples


plot_name <- "PCA1"
file_name <- paste0(root, "_", plot_name, ".pdf")

pdf(file_name)
plot(s$u[,1], s$u[,2], pch=19, cex=2, xlab=paste0("PC1, VarExp:", round(Var[1],3)), ylab=paste0("PC2, VarExp:", round(Var[2],3)),main="PCA plot TMM Normalization",col=col.cell, ylim=c((min(s$u[,2])*1.2),(max(s$u[,2])*1.2)),xlim=c((min(s$u[,1])*1.2),(max(s$u[,1])*1.2)))
text(s$u[,1], s$u[,2], labels=colnames(CPM),adj=c(0,2),cex=0.5)
legend("bottomright",legend = unique(metadata$groupA),col = unique(col.cell), pch=19,cex=1)
dev.off()

plot(s$u[,1], s$u[,2], pch=19, cex=2, xlab=paste0("PC1, VarExp:", round(Var[1],3)), ylab=paste0("PC2, VarExp:", round(Var[2],3)),main="PCA plot TMM Normalization",col=col.cell, ylim=c((min(s$u[,2])*1.2),(max(s$u[,2])*1.2)),xlim=c((min(s$u[,1])*1.2),(max(s$u[,1])*1.2)))
text(s$u[,1], s$u[,2], labels=colnames(CPM),adj=c(0,2),cex=0.5)
legend("bottomleft",legend = unique(metadata$groupA),col = unique(col.cell), pch=19,cex=1)






reads.cpm <- cpm(dge_norm)
unfiltered.results <- data.frame(id=dge_norm$genes, reads.cpm, dge.lrt$table)

filter <- apply(X=reads.cpm, MARGIN=1,
                FUN=function(data) {data[order(rank(data), decreasing=TRUE)[2]]})


lowerQuantile <- mean(filter == 0)
if (lowerQuantile < .95) upperQuantile <- .95 else upperQuantile <- 1
theta <- seq(lowerQuantile, upperQuantile, length=50)

filtPadj <- filtered_p(filter=filter, test=unfiltered.results$PValue,
                       theta=theta, method="BH")

min.fdr <- 0.05
numRej <- colSums(filtPadj < min.fdr, na.rm = TRUE)

filter.quantiles <- quantile(filter, probs=theta)

lo.fit.theta <- lowess(numRej ~ theta, f=1/5)

if (max(numRej) <= 10) {
  j <- 1
} else {
  residual <- if (all(numRej==0)) {
    0
  } else {
    numRej[numRej > 0] - lo.fit.theta$y[numRej > 0]
  }
  thresh <- max(lo.fit.theta$y) - sqrt(mean(residual^2))
  j <- if (any(numRej > thresh)) {
    which(numRej > thresh)[1]
  } else {
    1
  }
}

results <- unfiltered.results
results$FDR <- filtPadj[, j, drop=TRUE]

results$de <- sign(results$logFC)*(results$FDR < 0.05)
results$sig <- abs(results$de)

results[is.na(results$sig),]$sig <- 0
