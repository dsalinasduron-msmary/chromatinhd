# Make a small test dataset: a few clusters (cells), the most variable genes, the
# ATAC peaks near those genes, and a matching fragments file.
# Usage: Rscript make-test-subset.R <in.rds> <in.fragments.tsv.gz> <out.prefix> <n.genes> <flank.bp> <cluster_reordered> [<cluster_reordered> ...]
#   Clusters are the `cluster_reordered` numbers (see `cell_type_reordered_label`).
#   Genes: top <n.genes> of the SCT variable genes (scale.data) on chromosomes Gm01-Gm20,
#     ranked by SCT residual variance averaged over the per-library SCT models.
#   Peaks and fragments: kept within <flank.bp> of each gene (gene body +/- flank).
#   RNA is log-normalised over ALL genes before genes are dropped and stored as the
#     RNA data layer, so step 0 (lognorm) doesn't renormalise over the subset.
# Writes <out.prefix>.rds and <out.prefix>_fragments.tsv.gz (+ .tbi). Indexes the input
# fragments file (writes <in.fragments.tsv.gz>.tbi next to it) if it has no index.
suppressPackageStartupMessages({
	library(Seurat)
	library(Signac)
	library(GenomicRanges)
	library(Rsamtools)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 6) {
	stop("usage: Rscript make-test-subset.R <in.rds> <in.fragments.tsv.gz> <out.prefix> <n.genes> <flank.bp> <cluster_reordered> [...]")
}
in.rds <- args[1]
in.fragments <- args[2]
out.prefix <- args[3]
n.genes <- as.integer(args[4])
flank.bp <- as.integer(args[5])
clusters.oi <- as.integer(args[-(1:5)])
out.rds <- paste0(out.prefix, ".rds")
out.fragments <- paste0(out.prefix, "_fragments.tsv.gz")

obj <- readRDS(in.rds)

# cells
cells.oi <- colnames(obj)[obj$cluster_reordered %in% clusters.oi]

# genes: SCT variable genes ranked by mean residual variance over the SCT models
sct.models <- slot(obj[["SCT"]], "SCTModel.list")
res.var <- sapply(sct.models, function(m) {
	fa <- SCTResults(m, slot="feature.attributes")
	fa[rownames(obj[["SCT"]]), "residual_variance"]
})
rownames(res.var) <- rownames(obj[["SCT"]])
res.var <- rowMeans(res.var, na.rm=TRUE)
annotation <- Annotation(obj[["ATAC"]])
gene.annotation <- annotation[annotation$type == "gene"]
genes.main <- gene.annotation$gene_id[grepl("^Gm[0-9]+$", as.character(seqnames(gene.annotation)))]
genes.var <- rownames(GetAssayData(obj, assay="SCT", layer="scale.data"))
genes.var <- intersect(genes.var, genes.main)
genes.oi <- head(genes.var[order(res.var[genes.var], decreasing=TRUE)], n.genes)

# regions: gene bodies +/- flank; peaks overlapping them
gene.ranges <- gene.annotation[match(genes.oi, gene.annotation$gene_id)]
regions <- reduce(trim(resize(gene.ranges, width=width(gene.ranges) + 2*flank.bp, fix="center")), ignore.strand=TRUE)
start(regions) <- pmax(start(regions), 1)
peaks <- rownames(obj[["ATAC"]]) # granges() is in the same order
peaks.oi <- peaks[overlapsAny(granges(obj[["ATAC"]]), regions, ignore.strand=TRUE)]

# RNA: log-normalise over all genes first (library size = all genes, not the subset)
obj <- subset(obj, cells=cells.oi)
obj[["RNA"]] <- NormalizeData(JoinLayers(obj[["RNA"]]), verbose=FALSE)
obj.small <- subset(obj, features=c(genes.oi, peaks.oi))
obj.small@meta.data <- droplevels(obj.small@meta.data)

print(table(obj.small$cell_type_reordered_label, obj.small$libraries))
for (a in Assays(obj.small)) cat(a, ":", nrow(obj.small[[a]]), "features\n")
saveRDS(obj.small, out.rds)

# fragments: test cells within the regions
if (!file.exists(paste0(in.fragments, ".tbi"))) indexTabix(in.fragments, format="bed")
fragments.oi <- unique(unlist(scanTabix(TabixFile(in.fragments), param=regions), use.names=FALSE))
fragments.barcode <- vapply(strsplit(fragments.oi, "\t", fixed=TRUE), `[`, "", 4)
fragments.oi <- fragments.oi[fragments.barcode %in% cells.oi]
cat("fragments kept:", length(fragments.oi), "\n")
tmp.fragments <- paste0(out.prefix, "_fragments.tsv")
writeLines(fragments.oi, tmp.fragments)
bgzip(tmp.fragments, dest=out.fragments, overwrite=TRUE)
invisible(file.remove(tmp.fragments))
indexTabix(out.fragments, format="bed")
