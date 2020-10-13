# Some required R packages
library(Seurat)
library(edgeR)

# Code adapted from https://osca.bioconductor.org/multi-sample-comparisons.html#differential-expression-between-conditions

DEtest <- function(scObject, thre=100, outname){

# Construct pseudo bulk count matrix from Seurat Object
countlist <- lapply(names(table(scObject$SampleName)), 

				function(x){

					index <- which(scObject$SampleName == x)

					tmp <- scObject[["RNA"]]@counts[,index]

					if ( length(index) == 1){

						return (tmp)

						}else{

							return ( apply(tmp,1,sum) )

							}

					})

countmatrix <- do.call(cbind,countlist)

colnames(countmatrix) <- names(table(scObject$SampleName))

# Filter low expression genes
keep <- which(rowSums(countmatrix) > thre)

countmatrix <- countmatrix[keep,]

extra.info <- scObject@meta.data[match(colnames(countmatrix), scObject$SampleName),]

extra.info$Condition <- factor(extra.info$Condition, levels=c("control","treat"))

# Creating up a DGEList object for use in edgeR
y <- DGEList(countmatrix, samples=extra.info)

# Correct for composition biases by computing normalization factors
y <- calcNormFactors(y)

# Design matrix
design <- model.matrix(~factor(Condition), y$samples)

design

# Estimate the NB dispersion
y <- estimateDisp(y, design)

# Estimate the QL dispersion
fit <- glmQLFit(y, design, robust=TRUE)

# Test for differences in expression
res <- glmQLFTest(fit, coef=ncol(design))

saveRDS(res, file=outname)

}
