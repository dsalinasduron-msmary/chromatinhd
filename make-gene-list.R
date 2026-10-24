# Make a step-3 gene list (GENE_FILE in step-3.sh): the SCT variable genes (SCT scale.data)
# that step 2 sets up, ranked by SCT residual variance averaged over the per-library SCT
# models (highest first), as make-test-subset.R ranks them.
# Usage: Rscript make-gene-list.R <in.rds> <transcripts.csv> <out.genes.tsv>
#   transcripts.csv: step 1 output (TRANSCRIPTS in paths.sh); genes without a row there
#     (e.g. scaffold genes) are not set up by step 2 and are left out.
#   out.genes.tsv: a GENE header, then GLYMA_XXGXXXXXX names (the names step 2 uses).
# genes_sct-variable-2994.tsv was made with this from the cluster-reordered .rds.
suppressPackageStartupMessages({
	library(Seurat)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
	stop("usage: Rscript make-gene-list.R <in.rds> <transcripts.csv> <out.genes.tsv>")
}
in.rds <- args[1]
in.transcripts <- args[2]
out.genes <- args[3]

obj <- readRDS(in.rds)

# SCT variable genes ranked by mean residual variance over the SCT models
sct.models <- slot(obj[["SCT"]], "SCTModel.list")
res.var <- sapply(sct.models, function(m) {
	fa <- SCTResults(m, slot="feature.attributes")
	fa[rownames(obj[["SCT"]]), "residual_variance"]
})
rownames(res.var) <- rownames(obj[["SCT"]])
res.var <- rowMeans(res.var, na.rm=TRUE)
genes.var <- rownames(GetAssayData(obj, assay="SCT", layer="scale.data"))
genes.var <- genes.var[order(res.var[genes.var], decreasing=TRUE)]

# names as step 2 makes them: Glyma.09G186700.Wm82.a4.v1 --> GLYMA_09G186700
symbols <- toupper(gsub(".", "_", sub(".Wm82.a4.v1", "", genes.var, fixed=TRUE), fixed=TRUE))
symbols <- symbols[symbols %in% read.csv(in.transcripts)$symbol]

write.table(data.frame(GENE=symbols), out.genes, sep="\t", quote=FALSE, row.names=FALSE)
cat(length(genes.var), "SCT variable genes,", length(symbols), "with a TSS in", in.transcripts, "->", out.genes, "\n")
