# based off of https://chromatinhd.org/quickstart/1_data/#transcriptomics
# Step 2: build the ChromatinHD data (regions, fragments, transcriptome, folds).
# Usage: python step-2_training-setup.py <multiome.h5mu> <transcripts.csv> <fragments.tsv.gz> <out dir> <window bp>
#   multiome.h5mu:    paired RNA/ATAC from step 0 (expression to model in mod/RNA .X)
#   transcripts.csv:  one TSS per gene from step 1
#                     (chrom,start,end,transcript,strand,ensembl_gene_id,tss,symbol,transcript_biotype)
#   fragments.tsv.gz: ATAC fragments, bgzipped and tabix-indexed (.tbi)
#   out dir:          gets regions/, fragments/, transcriptome/, folds/5x1/ (input to steps 3-5;
#                     can use a lot of space) and genes.tsv (column GENE: every gene set up,
#                     the default gene list for step 3)
#   window bp:        region = TSS +/- window (e.g. 10000)
import sys
import pathlib

import chromatinhd as chd
import mudata
import pandas as pd

if len(sys.argv) != 6:
    sys.exit("usage: python step-2_training-setup.py <multiome.h5mu> <transcripts.csv> <fragments.tsv.gz> <out dir> <window bp>")
multiomic_file, transcripts_file, fragments_file, out_dir = sys.argv[1:5]
window_size = int(sys.argv[5])
out_dir = pathlib.Path(out_dir)
out_dir.mkdir(parents=True, exist_ok=True)

#########################
# BEGIN DATA PROCESSING #
#########################

mdata = mudata.read_h5mu(multiomic_file)
adata_rna = mdata["RNA"]
adata_atac = mdata["ATAC"]

# prep RNA data; change gene names and subset if necessary
genes_of_interest =  [ g for g in adata_rna.var.index ] # subset genes if needed here
rna_train = adata_rna[:,genes_of_interest] # grab only selected subset
# change names from  Glyma.09G186700.Wm82.a4.v1 --> GLYMA_09G186700 to maximize compatibility
ensembl_idx = [t.replace(".Wm82.a4.v1","").upper().replace(".","_") for t in genes_of_interest ]
rna_train.var_names = ensembl_idx

# generate transcripts database to cross-reference with fragments and rna expression data
transcripts = pd.read_csv(transcripts_file)
# remove scaffold genes
transcripts.drop(transcripts.index[ [ "scaffold" in ch for ch in transcripts["chrom"] ] ], inplace=True)

# can only train models for which we have transcript annotations
# subset training data to only those genes in transcripts
take = rna_train.var.index.intersection(transcripts.symbol)
transcripts.index = transcripts.symbol

transcripts_little = transcripts.loc[ take ]
rna_train_little = rna_train[:, take ]

rna_train = rna_train_little  # ignore the larger dataset from now on
transcripts = transcripts_little # ignore the larger dataset from now on

# index names become dimension names in ChromatinHD results (e.g. get_prediction
# fails on unnamed ones); match the transcriptome's "cell" / "gene"
transcripts.index.name = "gene"
adata_atac.obs.index.name = "cell"

regions = chd.data.Regions.from_transcripts(
    transcripts,
    path=out_dir / "regions",
    window=[-window_size, window_size]
)

fragments = chd.data.Fragments.from_fragments_tsv(
    fragments_file,
    regions,
    obs=adata_atac.obs,
    path=out_dir / "fragments",
)
fragments.create_regionxcell_indptr()

folds = chd.data.folds.Folds(out_dir / "folds" / "5x1").sample_cells(fragments, 5, 1)

# the transcriptome stores .X as layer "X" (lognorm or SCT, from step 0) and the raw
# counts as layer "counts"; step 3 trains on "X"
transcriptome = chd.data.Transcriptome.from_adata(rna_train, path=out_dir / "transcriptome")

# every gene set up here, in step 3's gene-list format
pd.DataFrame({"GENE": take}).to_csv(out_dir / "genes.tsv", sep="\t", index=False)

print(f"{len(take)} genes, {fragments.n_cells} cells, {fragments.coordinates.shape[0]} fragments -> {out_dir}")
