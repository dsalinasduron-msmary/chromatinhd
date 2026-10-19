# based off of https://chromatinhd.org/quickstart/1_data/#transcriptomics
# Step 3: train ChromatinHD prediction models (one per gene and fold).
# Usage: python step-3_train-models.py <data dir> <genes.tsv> <models dir>
#   data dir:   output of step 2 (regions/, fragments/, transcriptome/, folds/5x1/)
#   genes.tsv:  the genes to train models for, column GENE (GLYMA_XXGXXXXXX names)
#   models dir: models are added here (not reset, so step-3.sh can train in batches;
#               genes/folds already trained are skipped); can use a lot of space
import sys

import chromatinhd as chd
import chd_patches  # fixes for chromatinhd 0.4.3 (torch_scatter); see chd_patches.py
import pandas as pd

if len(sys.argv) != 4:
    sys.exit("usage: python step-3_train-models.py <data dir> <genes.tsv> <models dir>")
data_dir, genes_file, models_dir = sys.argv[1:]

fragments = chd.data.Fragments(data_dir + "/fragments")
transcriptome = chd.data.Transcriptome(data_dir + "/transcriptome")
folds = chd.data.folds.Folds(data_dir + "/folds/5x1")

# regions of interest = genes we want to train models for
# however, can only train models for genes in transcriptome data
genes = pd.read_csv(genes_file, sep="\t")["GENE"]
regions_oi = [g for g in genes if g in set(transcriptome.var.index)]
missing = len(genes) - len(regions_oi)
if missing:
    print(f"{missing} of {len(genes)} genes are not in the transcriptome and are skipped")

# Note: reset is set to false to enable batching, we batch to avoid "too many files open" error
models = chd.models.pred.model.multiscale.Models(models_dir, reset=False)
# train on the transcriptome's "X" layer (lognorm or SCT from step 0), not raw "counts";
# without this the model uses whichever layer happens to be listed first
models.model_params = dict(layer="X")

models.train_models(
    fragments=fragments,
    transcriptome=transcriptome,
    folds=folds,
    regions_oi=regions_oi
)
