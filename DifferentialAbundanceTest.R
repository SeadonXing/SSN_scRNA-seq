# Some required R packages
library(Seurat)
library(edgeR)

# Code adapted from https://osca.bioconductor.org/multi-sample-comparisons.html#differential-abundance

DAtest <- function(scObject, outname){

# Seurat Object
abundances <- table(scObject$FinalCelltype, scObject$SampleName) 

abundances <- unclass(abundances) 

head(abundances)

extra.info <- scObject@meta.data[match(colnames(abundances), scObject$SampleName),]

extra.info$Condition <- factor(extra.info$Condition, levels=c("control","treat"))

# Creating up a DGEList object for use in edgeR
y.ab <- DGEList(abundances, samples=extra.info)

y.ab

# Design matrix
design <- model.matrix(~factor(Condition), y.ab$samples)

# Estimate the NB dispersion
y.ab <- estimateDisp(y.ab, design, trend="none")

summary(y.ab$common.dispersion)

plotBCV(y.ab, cex=1)

# Estimate the QL dispersion
fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)

summary(fit.ab$var.prior)

summary(fit.ab$df.prior)

plotQLDisp(fit.ab, cex=1)

# Test for differences in abundance
res <- glmQLFTest(fit.ab, coef=ncol(design))

saveRDS(res, file=outname)

}
