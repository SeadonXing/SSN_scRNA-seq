# Some required R packages
library(Matrix)
library(Biobase)
library(Seurat)

CelltypeInteraction <- function(scObject, ligandReceptorFile, Eathreshold=0.5, Epthreshold=0.1, 
					out.detail="interaction.detail", out.stat="interaction.stat", Quanthreshold=3/4, ifOutFiles=TRUE){

	LRpairs <- ligandReceptorFile[,"Pair.Name"]

	ligands <- ligandReceptorFile[,"Ligand.ApprovedSymbol"]

	receptors <- ligandReceptorFile[,"Receptor.ApprovedSymbol"]

	# Seurat Object
	res <- lapply( levels(scObject$FinalCelltype), function(x){

		cell.used <- rownames(scObject@meta.data)[which(scObject$FinalCelltype==x)]

		rawE <- as.matrix(scObject[["RNA"]]@data)[, cell.used]

		geneEa <- apply(rawE,1,function(y){ log2(mean(exp(y)-1)+1) })

		geneEp <- apply(rawE,1,function(z){ length(which(z>0))/length(z) })

		tmp <- data.frame( geneName=rownames(rawE), Cluster=x, geneEa=geneEa, geneEp=geneEp, LorR="others", stringsAsFactors=FALSE )

		tmp[ tmp$geneName %in% ligands,"LorR" ] <- "Ligand"

		tmp[ tmp$geneName %in% receptors,"LorR" ] <- "Receptor"

		return (tmp)

	})

	inputdata <- do.call(rbind, res)

	inputdata <- inputdata[ (inputdata$geneEa >= Eathreshold) & (inputdata$geneEp >= Epthreshold), ]

	types <- unique(inputdata$Cluster)

	out <- lapply( 1:length(types), function(x){

			out.tmp <- lapply( x:length(types), function(y){

			x <- types[x]

			y <- types[y]

			print (paste0("Running ",x," and ",y," now !"))

			Exp.X <- inputdata[ inputdata$Cluster==x, ]

			Exp.Y <- inputdata[ inputdata$Cluster==y, ]

			Ligand.X <- Exp.X[ Exp.X$LorR=="Ligand", "geneName"]

			Receptor.X <- Exp.X[ Exp.X$LorR=="Receptor", "geneName"]

			Ligand.Y <- Exp.Y[ Exp.Y$LorR=="Ligand", "geneName"]

			Receptor.Y <- Exp.Y[ Exp.Y$LorR=="Receptor", "geneName"]

			#
			Lx_Ry <- lapply(Ligand.X,function(lx){ 

						tmp <- lapply(Receptor.Y,function(ry){ 

						 c( Ligand2Receptor=paste0(x," -> ",y),

							PairName_L_R=paste0(lx,"_",ry),

							LigandEa=Exp.X[ Exp.X$geneName==lx, "geneEa"],

							LigandEp=Exp.X[ Exp.X$geneName==lx, "geneEp"],

							ReceptorEa=Exp.Y[ Exp.Y$geneName==ry, "geneEa"],

							ReceptorEp=Exp.Y[ Exp.Y$geneName==ry, "geneEp"] ) }) 

						res <- as.data.frame(do.call(rbind,tmp))

					})

			Lx_Ry_tmp <- as.data.frame(do.call(rbind,Lx_Ry))

			Lx_Ry_real <- Lx_Ry_tmp[ Lx_Ry_tmp$PairName %in% LRpairs, ]

			#
			Ly_Rx <- lapply(Ligand.Y,function(ly){ 

						tmp <- lapply(Receptor.X,function(rx){ 

						 c( Ligand2Receptor=paste0(y," -> ",x),

							PairName_L_R=paste0(ly,"_",rx),

							LigandEa=Exp.Y[ Exp.Y$geneName==ly, "geneEa"],

							LigandEp=Exp.Y[ Exp.Y$geneName==ly, "geneEp"],

							ReceptorEa=Exp.X[ Exp.X$geneName==rx, "geneEa"],

							ReceptorEp=Exp.X[ Exp.X$geneName==rx, "geneEp"] ) }) 

						res <- as.data.frame(do.call(rbind,tmp))

					})

			Ly_Rx_tmp <- as.data.frame(do.call(rbind,Ly_Rx))

			Ly_Rx_real <- Ly_Rx_tmp[ Ly_Rx_tmp$PairName %in% LRpairs, ]

			return ( rbind(Lx_Ry_real, Ly_Rx_real))

		})

		return( as.data.frame(do.call(rbind,out.tmp)) )

	})

	out.final <- as.data.frame(do.call(rbind,out))

	if ( ifOutFiles ){

		write.table( out.final, paste0(out.detail,".txt"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")}

	#
	LigandNode <- unlist(lapply(strsplit(names(table(out.final$Ligand2Receptor))," -> "),function(m){m[1]}))

	ReceptorNode <- unlist(lapply(strsplit(names(table(out.final$Ligand2Receptor))," -> "),function(n){n[2]}))

	stat.raw <- data.frame( "LigandNode"= LigandNode,

							"ReceptorNode"= ReceptorNode,

							"LRnum"=as.numeric(table(out.final$Ligand2Receptor)),

							"LigangAttri"=as.numeric(table(inputdata$Cluster)[LigandNode]),

							"ReceptorAttri"=as.numeric(table(inputdata$Cluster)[ReceptorNode]))

	if ( ifOutFiles ){

	write.table( stat.raw, paste0(out.stat,".raw.txt"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")}

	keepNum <- quantile(stat.raw$LRnum,Quanthreshold)

	stat.filter <- stat.raw[ stat.raw$LRnum >= keepNum, ]

	if ( ifOutFiles ){

	write.table( stat.filter, paste0(out.stat,".filter.txt"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")}

	return ( list(out.final=out.final, stat.raw=stat.raw, stat.filter=stat.filter) )

}
