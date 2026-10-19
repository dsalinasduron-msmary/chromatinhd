# Convert the integrated Seurat multiome object to .h5mu for ChromatinHD (step 2).
# Usage: Rscript step-0_convert-rds-to-h5mu.R <in.rds> <out.h5mu> <lognorm|sct>
#   lognorm: RNA modality = raw RNA counts, log-normalised here (NormalizeData),
#            or the RNA data layer if the object already has one (test subset)
#   sct:     RNA modality = the object's SCTransform-normalised data
# Either way the expression ChromatinHD uses ends up in mod/RNA .X; raw counts
# are kept in layers['counts'].
suppressPackageStartupMessages({
	library(Seurat)
	library(Signac)
	library(MuDataSeurat)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3 || !(args[3] %in% c("lognorm","sct"))) {
	stop("usage: Rscript step-0_convert-rds-to-h5mu.R <in.rds> <out.h5mu> <lognorm|sct>")
}
in.rds <- args[1]
out.h5mu <- args[2]
method <- args[3]

integrated.multiome.rds <- readRDS(in.rds)

# MuDataSeurat cannot write Seurat v5 (Assay5) assays, so rebuild the RNA
# modality as a v3 Assay holding counts + log-normalised data (WriteH5MU
# then writes data to .X and counts to layers['counts']).
if (method == "lognorm") {
	rna <- JoinLayers(integrated.multiome.rds[["RNA"]])
	rna.assay <- CreateAssayObject(counts=LayerData(rna, layer="counts"))
	# use an existing data layer (the gene-subset test data, normalised over all genes
	# by make-test-subset.R); otherwise normalise here
	if ("data" %in% Layers(rna)) {
		rna.assay <- SetAssayData(rna.assay, layer="data", new.data=LayerData(rna, layer="data"))
	} else {
		rna.assay <- NormalizeData(rna.assay, verbose=FALSE)
	}
} else {
	sct <- integrated.multiome.rds[["SCT"]]
	rna.assay <- CreateAssayObject(counts=LayerData(sct, layer="counts"))
	rna.assay <- SetAssayData(rna.assay, layer="data", new.data=LayerData(sct, layer="data"))
}
DefaultAssay(integrated.multiome.rds) <- "ATAC"
integrated.multiome.rds[["SCT"]] <- NULL # its scale.data (variable genes only) makes an invalid modality
integrated.multiome.rds[["RNA"]] <- rna.assay
DefaultAssay(integrated.multiome.rds) <- "RNA"
# WNN graphs record no assay, which crashes WriteH5MU; ChromatinHD doesn't use them
integrated.multiome.rds@graphs <- list()

if (file.exists(out.h5mu)) invisible(file.remove(out.h5mu))
tryCatch(
	WriteH5MU(integrated.multiome.rds, out.h5mu),
	error=function(e) {
		if (file.exists(out.h5mu)) invisible(file.remove(out.h5mu)) # don't leave a partial file behind
		stop(e)
	}
)
