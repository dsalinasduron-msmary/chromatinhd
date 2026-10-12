# based off of https://chromatinhd.org/quickstart/1_data/#transcriptomics
import torch
import anndata as ad
import chromatinhd as chd
import mudata
import pandas as pd

################################
# FILES PATHS  & OTHER OPTIONS #
# CHANGE AS NEEDED             #
################################

pt = "/falafel/cornbert" # path to put code artifacts (including output models, can consume a lot of space)

multiomic_file = "/falafel/cornbert/integrated.multiome.rds.h5mu" # paired atac/rna data

# transcription start site spec: chrom,start,end,transcript,strand,ensembl_gene_id,tss,symbol,transcript_biotype
transcripts_file = "/falafel/cornbert/transcripts.csv" 

fragments_file = "/falafel/pgrp_singlecell/soybean_nodule_multiome-redo/atac_fragments.tsv.gz"

model_label = "test" # will be added as suffix to artifacts
window_size = 10_000

# set to true generate index if fragments does not have one [error will tip you off if you have to do this]
if False:
    import pysam
    pysam.tabix_index(fragments_file, preset = "bed")

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

regions = chd.data.Regions.from_transcripts(
    transcripts,
    path=(pt + "/" + "regions_" + model_label),
    window=[-window_size, window_size]
)

fragments = chd.data.Fragments.from_fragments_tsv(
    fragments_file,
    regions,
    obs=adata_atac.obs,
    path= (pt + "/" + "fragments_" + model_label),
)
fragments.create_regionxcell_indptr()

folds = chd.data.folds.Folds(pt + "/" + "folds_" + model_label  + "/5x1").sample_cells(fragments, 5, 1)

transcriptome = chd.data.Transcriptome.from_adata(rna_train, path=(pt + "/" + "transcriptome_" + model_label) )
