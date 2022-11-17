#
# Builds Seurat object used for downstream analyses
#
library(Matrix)
library(Seurat)
suppressMessages(library(tidyverse))

cfg <- snakemake@config

# load data, barcodes and feature metadata
barcodes <- read_tsv(cfg$input_data$cell_ids, col_names = c("barcode"), show_col_types = FALSE)

features <- read_tsv(cfg$input_data$genes,
                     col_names = c('ensgene', 'symbol', 'type'),
                     show_col_types = FALSE)

mat <- mat <- readMM(cfg$input_data$counts)
rownames(mat) <- features$ensgene
colnames(mat) <- barcodes$barcode

# construct seurat object
seurat_object <- CreateSeuratObject(counts = mat, project = "scRNA-seq", assay = "RNA")

# save as rds
saveRDS(seurat_object, snakemake@output[[1]])
