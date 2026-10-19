# Step 1a: candidate transcripts (TSSs) for ChromatinHD from the gene annotation
# embedded in the Seurat/Signac object (Annotation(obj[["ATAC"]]), Phytozome v13).
# Writes every transcript of every gene in the RNA assay; step-1_select-tss.py then
# picks one TSS per gene from the fragments and writes transcripts.csv for step 2.
# Usage: Rscript step-1_annotation-to-transcripts.R <in.rds> <out.candidates.csv>
# Coordinates are 1-based (as in the annotation); strand is 1/-1 (ChromatinHD);
# tss = start on the + strand, end on the - strand.
suppressPackageStartupMessages({
	library(Seurat)
	library(Signac)
	library(GenomicRanges)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
	stop("usage: Rscript step-1_annotation-to-transcripts.R <in.rds> <out.candidates.csv>")
}
in.rds <- args[1]
out.csv <- args[2]

integrated.multiome.rds <- readRDS(in.rds)
annotation <- Annotation(integrated.multiome.rds[["ATAC"]])
genes.rna <- rownames(integrated.multiome.rds[["RNA"]])

transcripts <- annotation[annotation$type == "transcript" & annotation$gene_id %in% genes.rna]
stopifnot(all(as.character(strand(transcripts)) %in% c("+", "-")))
transcripts.df <- data.frame(
	chrom=as.character(seqnames(transcripts)),
	start=start(transcripts),
	end=end(transcripts),
	strand=ifelse(as.character(strand(transcripts)) == "+", 1, -1),
	transcript=transcripts$transcript_id,
	ensembl_gene_id=transcripts$gene_id,
	gene_biotype=transcripts$gene_biotype
)
transcripts.df$tss <- ifelse(transcripts.df$strand == 1, transcripts.df$start, transcripts.df$end)
transcripts.df <- unique(transcripts.df)

genes.missing <- setdiff(genes.rna, transcripts.df$ensembl_gene_id)
cat("RNA genes:", length(genes.rna), " with transcripts:", length(unique(transcripts.df$ensembl_gene_id)),
	" candidate transcripts:", nrow(transcripts.df), "\n")
if (length(genes.missing)) cat("RNA genes without a transcript in the annotation:", length(genes.missing), "\n")

write.csv(transcripts.df, out.csv, row.names=FALSE, quote=FALSE)
