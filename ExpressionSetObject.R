# Some required R packages
library(Seurat)
library(dplyr)
library(Matrix)
library(Biobase)

ExpressionSetObject <- function(scObject, ifRawCounts=TRUE, annotation.name, outname){

# Check Seurat Object
DefaultAssay(object = scObject)

levels(Idents(scObject))

head(scObject@meta.data)

# assayData
if (ifRawCounts){

assayData <- scObject[["RNA"]]@counts # raw counts

}else{

assayData <- scObject[["RNA"]]@data # normalized data

}

dim(assayData)

assayData[1:3,1:3]

# phenoData
cell.anno <- scObject@meta.data

cell.metaname <- data.frame(labelDescription=colnames(cell.anno),

				row.names=colnames(cell.anno));

phenoData <- new("AnnotatedDataFrame",data=cell.anno,varMetadata=cell.metaname)

phenoData

# featureData
gene.anno <- data.frame( gene_short_name=rownames(assayData),

			Database="10XEnsmble",row.names=rownames(assayData))

gene.metaname <- data.frame(labelDescription=

			c('Symbols of Genes', 'Gene symbol original database'),

			row.names=c('gene_short_name', 'Database'));

featureData <-new("AnnotatedDataFrame",data=gene.anno,varMetadata=gene.metaname)

featureData

# Check
print (identical( colnames(assayData), rownames(cell.anno)))

print (identical( rownames(assayData), rownames(gene.anno)))

# ExpressionSet Object
ExpressionSetObject <- ExpressionSet(assayData= as.matrix(assayData),

						phenoData=phenoData,

						featureData=featureData,

						annotation=annotation.name)

ExpressionSetObject

# Save RDS
saveRDS(ExpressionSetObject, file = outname);

}
