# Make a small test dataset by keeping only a few clusters of the full object.
# Usage: Rscript make-test-subset.R <in.rds> <out.rds> <cluster_reordered> [<cluster_reordered> ...]
# Clusters are the `cluster_reordered` numbers (see `cell_type_reordered_label`).
suppressPackageStartupMessages({
	library(Seurat)
	library(Signac)
})

args <- commandArgs(trailingOnly=TRUE)
in.rds <- args[1]
out.rds <- args[2]
clusters.oi <- as.integer(args[-(1:2)])

obj <- readRDS(in.rds)
cells.oi <- colnames(obj)[obj$cluster_reordered %in% clusters.oi]
obj.small <- subset(obj, cells=cells.oi)
obj.small@meta.data <- droplevels(obj.small@meta.data)

print(table(obj.small$cell_type_reordered_label, obj.small$libraries))
saveRDS(obj.small, out.rds)
