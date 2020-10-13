# Required R packages
library("infercnv")

# Code adapted from : https://github.com/broadinstitute/inferCNV/wiki

# create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="inferCNV.counts.matrix",

                                    annotations_file="inferCNV.annotations.txt",

                                    delim="\t",

                                    gene_order_file="gencode_v19_gene_pos_rmMTXY.txt",

                                    ref_group_names=c("T cells","B cells","Plasma cells")) # names of immune cells

# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,

                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics

                             min_cells_per_gene = 3,

                             out_dir="output_dir", 

                             cluster_by_groups=T,

                             plot_steps=F,

                             mask_nonDE_genes = T,

                             k_obs_groups = 1,

                             hclust_method = "complete",

                             include.spike=F # used for final scaling to fit range (0,2) centered at 1.

)
