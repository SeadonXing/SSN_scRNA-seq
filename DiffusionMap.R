# Some required R packages
library(statmod)
library(destiny)
library(Matrix)
library(rgl)
library(Seurat)
library(dplyr)

DiffusionMap <- function(scObject, celltype.vector, variable.gene.num, outname){

# Check Seurat Object
DefaultAssay(object = scObject)

levels(Idents(scObject))

# Subset cells
scObject <- subset(x = scObject, 

			idents = celltype.vector, invert = FALSE)

# Get variable genes
scObject.markers <- FindAllMarkers(object = scObject, only.pos = TRUE, min.pct = 0.2, 

					logfc.threshold = 0.25, max.cells.per.ident=200, assay="RNA", slot = "data")

tmp <- as.data.frame(scObject.markers %>% group_by(cluster) %>% top_n(n = variable.gene.num, wt = avg_logFC))

vg <- tmp$gene

# Diffusion maps
runData <- as.matrix(scObject[["RNA"]]@data[ vg, ])

dm <- DiffusionMap(t(runData))

dif <- data.frame(dm@eigenvectors)

rownames(dif) = colnames(runData)

res <- cbind(dif, scObject@meta.data)

saveRDS(res, file=outname)

}
