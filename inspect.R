# Print a Seurat/Signac object's assays, layers, metadata, clusters, barcodes and fragment paths
# (e.g. to check an object on the production server before step 0).
# Usage: Rscript inspect.R <object.rds>
suppressPackageStartupMessages({
	library(Seurat)
	library(Signac)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {
	stop("usage: Rscript inspect.R <object.rds>")
}
in.rds <- args[1]

start.time <- Sys.time()
obj <- readRDS(in.rds)
cat("load time:", format(Sys.time() - start.time), " object size:", format(object.size(obj), units="GB"), "\n")
print(obj)

cat("\nAssays:\n")
for (assay in Assays(obj)) {
	layers <- tryCatch(paste(Layers(obj[[assay]]), collapse=","), error=function(e) "?")
	cat(assay, class(obj[[assay]])[1], "layers:", layers, " dim:", dim(obj[[assay]]), "\n")
}
cat("\nDefault assay:", DefaultAssay(obj), "\nReductions:", Reductions(obj), "\nIdents levels:", nlevels(Idents(obj)), "\n")

cat("\nmeta.data columns:\n")
str(obj@meta.data, list.len=60, vec.len=2)
print(table(Idents(obj)))

cat("\nBarcodes:\n")
print(head(colnames(obj)))
print(tail(colnames(obj), 3))

for (assay in Assays(obj)) {
	if (!inherits(obj[[assay]], "ChromatinAssay")) next
	cat("\nChromatinAssay", assay, "fragments objects:", length(Fragments(obj[[assay]])), "\n")
	for (fragment in Fragments(obj[[assay]])) {
		cat(GetFragmentData(fragment, "path"), " cells:", length(Cells(fragment)), "\n")
		print(head(GetFragmentData(fragment, "cells"), 3))
	}
	cat("links:", length(Links(obj[[assay]])), " annotation:", length(Annotation(obj[[assay]])), "\n")
	print(head(rownames(obj[[assay]]), 3))
}

cat("\nRNA features:\n")
print(head(rownames(obj[["RNA"]]), 3))
