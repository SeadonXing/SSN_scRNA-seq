# Some required R packages
library(Seurat)
library(dplyr)
library(Matrix)
library(GSVA)
library(GSEABase)
library(GSVAdata)
library(limma)
library(msigdbr)
library(Biobase)
library(pheatmap)

# Code adapted from : http://www.bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA.pdf
#					  http://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf

GSVAfunction <- function(ExpressionSetObject, celltype.vector, 
			Exp=0.1, Pct=0.01, GmtFile, min.size=5, max.size=500, outname){

# some checks for ExpressionSet Object
names(pData(ExpressionSetObject))
table(ExpressionSetObject$FinalCelltype) # The assigned cell types for each cell
table(ExpressionSetObject$Condition) # The group for each cell, control or treat

# select cells
selectExpressionSetObject <- ExpressionSetObject[ , ExpressionSetObject$FinalCelltype %in% celltype.vector ]

selectExpressionSetObject

# filter low expression genes
candidategenes <- rownames(exprs(selectExpressionSetObject))

length(candidategenes)

# log2 mean expression >= Exp

meanExp <- lapply( candidategenes, function(x) { 

	tmp <-as.data.frame(tapply(exprs(selectExpressionSetObject)[x,], 

		selectExpressionSetObject$FinalCelltype, 

		function(y){ log2(mean(exp(y)-1)+1) } ));

		colnames(tmp) <- x ;

	tmp 

} )

meanExpData <- as.data.frame(t(do.call(cbind,meanExp)))

saveGeneByMean <- rownames(meanExpData[ rowMeans(meanExpData )>= Exp ])

length(saveGeneByMean)

# mean express percent >= Pct
expPercent <- lapply( saveGeneByMean, function(x) { 

		tmp <-as.data.frame(tapply(exprs(selectExpressionSetObject)[x,], 

			selectExpressionSetObject$FinalCelltype, 

			function(y){ length(which(y>0))/length(y) } ));

			colnames(tmp) <- x ;

	tmp

} )

expPercentData <- as.data.frame(t(do.call(cbind,expPercent)))

saveGeneByExpPercent <- rownames(expPercentData[ rowMeans(expPercentData) >= Pct, ])

length(saveGeneByExpPercent)

# Filtered selectExpressionSetObject
selectExpressionSetObject <- selectExpressionSetObject[saveGeneByExpPercent,]

selectExpressionSetObject

# run GSVA
geneset <- getGmt(GmtFile)

GSVAresults <- gsva( selectExpressionSetObject , geneset , min.sz=min.size, max.sz=max.size,

			method="gsva", kcdf="Gaussian", mx.diff=TRUE, verbose=TRUE, parallel.sz=1)

GSVAresults

# save GSVA object
saveRDS(GSVAresults, file = paste0(outname,".GSVA.rds"))

# run Limma
condition <- factor( GSVAresults$Condition, levels=c("control","treat") )

# design matrix
design <- model.matrix(~ condition)

# fitting a linear model
fit <- lmFit(GSVAresults, design)

# apply empirical Bayes smoothing to the standard errors
fit <- eBayes(fit)

# save fit object
saveRDS(fit, file = paste0(outname,".fit.rds"))

}
