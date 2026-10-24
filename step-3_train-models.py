# based off of https://chromatinhd.org/quickstart/1_data/#transcriptomics
# Step 3: train ChromatinHD prediction models (one per gene and fold).
# Usage: python step-3_train-models.py <data dir> <genes.tsv> <models dir>
#   data dir:   output of step 2 (regions/, fragments/, transcriptome/, folds/5x1/)
#   genes.tsv:  the genes to train models for, column GENE (GLYMA_XXGXXXXXX names)
#   models dir: models are added here (not reset, so step-3.sh can train in batches;
#               genes/folds already trained are skipped); can use a lot of space
import math
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
model_params = dict(layer="X")

# Fragment encoder bin widths. chd's default spline encoder bins the region at 100, 200,
# 500, 1000, 2000 and 5000 bp = the 10 kb half-window divided by 2, 5, 10, 20, 50, 100
# (a 1-2-5 series). Each width must divide the window (TSS +/- WINDOW): otherwise the
# last, partial bin indexes past the encoder weights (CUDA device-side assert) or into
# the next width's weights (silently wrong), so the default widths only fit multiples
# of 10 kb. Use the same divisors of the actual half-window instead (3 kb: 30, 60, 150,
# 300, 600, 1500 bp): always six distinct widths that divide the window, in the same
# ratios and with the same number of encoder weights as the default; needs WINDOW to be
# a multiple of 100. The encoder is saved with each model, so steps 4-5 use the same widths.
default_binwidths = (100, 200, 500, 1000, 2000, 5000)
binwidth_divisors = (100, 50, 20, 10, 5, 2)
window = [int(w) for w in fragments.regions.window]
half_window = math.gcd(-window[0], window[1])  # widths must divide both ends of the window
if half_window % 100:
    sys.exit(f"window {window[0]}..{window[1]} bp: WINDOW must be a multiple of 100 bp "
             f"(encoder bin widths = WINDOW / {binwidth_divisors})")
binwidths = tuple(half_window // k for k in binwidth_divisors)
if binwidths != default_binwidths:
    print(f"WARNING: window {window[0]}..{window[1]} bp is not +/- 10 kb; using custom encoder "
          f"bin widths {binwidths} bp (WINDOW / {binwidth_divisors}) instead of chromatinhd's "
          f"default {default_binwidths}", file=sys.stderr)
    model_params["encoder_kwargs"] = dict(binwidths=binwidths)
else:
    print(f"encoder bin widths: {binwidths} bp (chromatinhd default)")
models.model_params = model_params

models.train_models(
    fragments=fragments,
    transcriptome=transcriptome,
    folds=folds,
    regions_oi=regions_oi
)
