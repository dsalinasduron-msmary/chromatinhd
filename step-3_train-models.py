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

gf =  "test-genes.tsv" #  the genes we want to train models for
pt = "/falafel/cornbert" # path to put code artifacts (including output models, can consume a lot of space)
model_label = "window-10kb"

fragments = chd.data.Fragments(pt + "/" + "fragments_" + model_label )
transcriptome = chd.data.Transcriptome(pt + "/" + "transcriptome_" + model_label)
folds = chd.data.folds.Folds(pt + "/" + "folds_" + model_label  + "/5x1")

# regions of interest = genes we want to train models for
# however, can only train models for genes in transcriptome data
take = set(transcriptome.adata.var.index)
region_oi_df = pd.read_csv(gf,sep="\t")
regions_oi = list( take.intersection(region_oi_df["GENE"]) )

# Note: reset is set to false to enable batching, we batch to avoid "too many files open" error
models = chd.models.pred.model.multiscale.Models(pt + "/" + "models_" + model_label, reset=False)

models.train_models(
    fragments=fragments,
    transcriptome=transcriptome,
    folds=folds,
    regions_oi=regions_oi
)
