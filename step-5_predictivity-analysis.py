# based off of https://chromatinhd.org/quickstart/3_pred/#which-regions-are-most-predictive
# Step 5: predictivity analysis - which parts of each gene's region the model relies on.
# For every model (gene x fold) windows of 50-500 bp are hidden one at a time and the drop
# in correlation between prediction and expression (deltacor) is measured.
# Usage: python step-5_predictivity-analysis.py <data dir> <models dir> <correlations.tsv>
#                                               <min cor> <min deltacor> <out dir>
#   data dir:         output of step 2 (regions/, fragments/, transcriptome/, folds/5x1/)
#   models dir:       output of step 3
#   correlations.tsv: per-gene table from step 4 (symbol, cor, cor_n_fragments, deltacor, n_folds)
#   min cor, min deltacor: genes with cor >= min cor and deltacor >= min deltacor (step-4
#                     means over folds) are analysed; genes with NaN values are not
#   out dir:          gets
#     selected_genes.tsv      the rows of correlations.tsv that passed
#     regionmultiwindow/      ChromatinHD results; genes already scored are skipped on reruns
#     predictivity.tsv        per gene, every 10 bp: symbol, position (bp from the TSS,
#                             + = downstream on the gene's strand), deltacor (< 0: hiding
#                             the window makes the prediction worse), lost (fragments
#                             hidden per cell), effect (> 0: hiding lowers the predicted
#                             expression, i.e. accessibility there goes with higher expression)
#     predictive_windows.tsv  stretches with deltacor < -0.001 (ChromatinHD default), split
#                             by effect_direction (+1 / -1); start/end relative to the TSS
#                             plus chrom/start_genome/end_genome
#     plots/<symbol>.png      lost fragments (pileup) and predictivity along the region
import sys
import pathlib

import chromatinhd as chd
import chd_patches  # fixes for chromatinhd 0.4.3 (torch_scatter, plot panels, fragment cache); see chd_patches.py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import polyptich

if len(sys.argv) != 7:
    sys.exit("usage: python step-5_predictivity-analysis.py <data dir> <models dir> <correlations.tsv> <min cor> <min deltacor> <out dir>")
data_dir, models_dir, correlations_file = sys.argv[1:4]
min_cor, min_deltacor = float(sys.argv[4]), float(sys.argv[5])
out_dir = pathlib.Path(sys.argv[6])
(out_dir / "plots").mkdir(parents=True, exist_ok=True)

fragments = chd.data.Fragments(data_dir + "/fragments")
transcriptome = chd.data.Transcriptome(data_dir + "/transcriptome")
folds = chd.data.folds.Folds(data_dir + "/folds/5x1")
models = chd.models.pred.model.multiscale.Models(models_dir, reset=False)

# genes to analyse: good predictions that beat counting fragments, with a model for every fold
gene_cors = pd.read_csv(correlations_file, sep="\t")
selected = gene_cors[(gene_cors["cor"] >= min_cor) & (gene_cors["deltacor"] >= min_deltacor)]
trained = [all(f"{g}_{i}" in models for i in range(len(folds))) for g in selected["symbol"]]
if not all(trained):
    print(f"{len(trained) - sum(trained)} selected genes lack a model for some fold and are skipped")
selected = selected[trained]
selected.to_csv(out_dir / "selected_genes.tsv", sep="\t", index=False)
genes = list(selected["symbol"])
print(f"{len(genes)} of {len(gene_cors)} genes with cor >= {min_cor} and deltacor >= {min_deltacor}")
if not genes:
    sys.exit("no genes selected; lower the thresholds")

# set up over all regions (results are sparse: only scored genes take space), so a
# rerun with other genes adds to it
censorer = chd.models.pred.interpret.MultiWindowCensorer(fragments.regions.window)
regionmultiwindow = chd.models.pred.interpret.RegionMultiWindow.create(
    path=out_dir / "regionmultiwindow",
    folds=folds,
    transcriptome=transcriptome,
    censorer=censorer,
    fragments=fragments,
    overwrite=False,
)
regionmultiwindow.score(
    models=models,
    regions=genes,
    folds=folds,
    fragments=fragments,
    transcriptome=transcriptome,
)
regionmultiwindow.interpolate(regions=genes)

predictivity = pd.concat(
    [regionmultiwindow.get_plotdata(symbol).reset_index().assign(symbol=symbol) for symbol in genes]
)
predictivity = predictivity[["symbol", "position", "deltacor", "lost", "effect"]]
predictivity.to_csv(out_dir / "predictivity.tsv", sep="\t", index=False)

windows = regionmultiwindow.extract_predictive_windows(region_id=genes).reset_index(drop=True)
windows = windows.rename(columns={windows.columns[2]: "symbol"})  # named after the region index ("gene")
windows["start"] = windows["start"].astype(int)
windows["end"] = windows["end"].astype(int)
windows = chd.data.regions.uncenter_multiple(windows.assign(region=windows["symbol"]), fragments.regions.coordinates)
windows = windows[["symbol", "start", "end", "effect_direction", "chrom", "start_genome", "end_genome"]]
windows.to_csv(out_dir / "predictive_windows.tsv", sep="\t", index=False)

for symbol in genes:
    fig = polyptich.grid.Figure(polyptich.grid.Grid(padding_height=0.05))
    fig.main.add_under(polyptich.grid.Title(symbol))
    fig.main.add_under(chd.models.pred.plot.Pileup.from_regionmultiwindow(regionmultiwindow, symbol, width=10))
    fig.main.add_under(chd.models.pred.plot.Predictivity.from_regionmultiwindow(regionmultiwindow, symbol, width=10))
    fig.plot()
    # matplotlib's savefig: polyptich's own imports IPython (to show it in a notebook)
    matplotlib.figure.Figure.savefig(fig, out_dir / "plots" / (symbol + ".png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

print(f"{len(genes)} genes, {len(windows)} predictive windows -> {out_dir}")
