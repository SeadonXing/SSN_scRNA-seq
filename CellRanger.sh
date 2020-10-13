# Code adapted from https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/2.2/using/count
#!/bin/bash

export PATH=$PATH:/Software/cellranger-2.2.0/cellranger-cs/2.2.0/bin

cellranger count --id=Cellranger_count \

                 --transcriptome=/CellRanger/refdata-cellranger-GRCh38-1.2.0 \

                 --fastqs=/CellRanger_fastq \

                 --expect-cells=5000 \

                 --localcores=32 \

                 --localmem=64
