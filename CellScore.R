# Some required R packages
library(Matrix)
library(Biobase)

# Code adapted from the 'AddModuleScore' function in Seurat R package (https://github.com/satijalab/seurat/blob/fe93b05745e55ec2f66e3f0b4c4196aad9f4d5a7/R/utilities.R).

CellScore <- function( inputdata, # the inputdata can be a "data.frame", "matrix" or "Sparse matrix"
			genes.list = NULL, genes.pool = NULL, n.bin = 25, ctrl.size = 100, 
			enrich.name = "Module_A", random.seed = 1){

# set seed
set.seed(seed = random.seed)

#transfor the inputdata to sparse matrix
if (inherits(x = inputdata, what = "data.frame")) {

    inputdata <- as.matrix(x = inputdata)

}

inputdata <- as(object = inputdata, Class = "RsparseMatrix")

# get genes in genes.list that exist in input data
genes.list <- lapply( X = genes.list, FUN = function(x) {

      return(intersect(x = x, y = rownames(x = inputdata)))})

cluster.length <- length(x = genes.list) # cluster.length=1

# in not specific background genes, use all the genes in data
if (is.null(x = genes.pool)) {

  genes.pool = rownames(x = inputdata)

}

# according to gene meanexp in all cells to split genes into 25 bins
data.avg <- Matrix::rowMeans(x = inputdata[genes.pool, ])

data.avg <- data.avg[order(data.avg)]

data.cut <- as.numeric(x = Hmisc::cut2(

  x = data.avg,

  m = round(x = length(x = data.avg) / n.bin)

))

names(x = data.cut) <- names(x = data.avg)

# for each input genes, get number of ctrl.size genes in the same bins, 
# finally get a control gene set with Num(input)*Ctrl.size
ctrl.use <- vector(mode = "list", length = cluster.length)

for (i in 1:cluster.length) {

  genes.use <- genes.list[[i]]

  for (j in 1:length(x = genes.use)) {

    ctrl.use[[i]] <- c(

      ctrl.use[[i]],

      names(x = sample(

        x = data.cut[which(x = data.cut == data.cut[genes.use[j]])],

        size = ctrl.size,

        replace = FALSE)))

}}

ctrl.use <- lapply(X = ctrl.use, FUN = unique)

# calculate ctrl gene set score
ctrl.scores <- matrix(

  data = numeric(length = 1L),

  nrow = length(x = ctrl.use),

  ncol = ncol(x = inputdata)

)

for (i in 1:length(ctrl.use)) {

  genes.use <- ctrl.use[[i]]

  ctrl.scores[i, ] <- Matrix::colMeans(x = inputdata[genes.use, ])

}

# calculate real module score
genes.scores <- matrix(

  data = numeric(length = 1L),

  nrow = cluster.length,

  ncol = ncol(x = inputdata)

)

for (i in 1:cluster.length) {

  genes.use <- genes.list[[i]]

  data.use <- inputdata[genes.use, , drop = FALSE]

  genes.scores[i, ] <- Matrix::colMeans(x = data.use)

}

# substract background ctrl genes score from real input module genes score
genes.scores.use <- genes.scores - ctrl.scores

if (length(cluster.length)==1){

    rownames(x = genes.scores.use) <- enrich.name

} else{

    rownames(x = genes.scores.use) <- paste0(enrich.name, 1:cluster.length)

}

genes.scores.use <- as.data.frame(x = t(x = genes.scores.use))

rownames(x = genes.scores.use) <- colnames(x = inputdata)

# return results
return(genes.scores.use)

}
