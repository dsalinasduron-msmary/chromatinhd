# based off of https://chromatinhd.org/quickstart/3_pred/#some-quality-checks
# Step 4: correlation analysis - how well each trained model predicts expression in its
# fold's test cells, compared with just counting fragments in the region.
# (Not predictivity analysis: that is step 5, which hides parts of the region.)
# Usage: python step-4_correlation-analysis.py <data dir> <models dir> <out dir>
#   data dir:   output of step 2 (regions/, fragments/, transcriptome/, folds/5x1/)
#   models dir: output of step 3; every trained model in it is scored
#   out dir:    gets
#     correlations_per-fold.tsv  one row per gene x fold: symbol, fold, cor (ChromatinHD
#                                prediction vs expression), cor_n_fragments (fragment
#                                count vs expression), deltacor = cor - cor_n_fragments.
#                                NaN where expression (or prediction / fragment count) is
#                                constant across the fold's test cells.
#     correlations.tsv           one row per gene: means over folds (NaN folds skipped),
#                                n_folds = folds with a valid cor. Filter this for the
#                                step-5 gene list (column "symbol").
#     correlation_plot.png       per-gene cor_n_fragments vs cor
import sys
import pathlib

import chromatinhd as chd
import chd_patches  # fixes for chromatinhd 0.4.3 (torch_scatter, fragment cache); see chd_patches.py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

if len(sys.argv) != 4:
    sys.exit("usage: python step-4_correlation-analysis.py <data dir> <models dir> <out dir>")
data_dir, models_dir, out_dir = sys.argv[1:]
out_dir = pathlib.Path(out_dir)
out_dir.mkdir(parents=True, exist_ok=True)

fragments = chd.data.Fragments(data_dir + "/fragments")
transcriptome = chd.data.Transcriptome(data_dir + "/transcriptome")
folds = chd.data.folds.Folds(data_dir + "/folds/5x1")
models = chd.models.pred.model.multiscale.Models(models_dir, reset=False)

# same computation as models.get_region_cors (chromatinhd 0.4.3), which drops the fold
# number; scores the full data directly (no per-gene filtering, which was slow)
def cor(x, y):
    with np.errstate(invalid="ignore", divide="ignore"):  # constant input -> NaN
        return np.corrcoef(x, y)[0, 1]

rows = []
for symbol in fragments.var.index:
    for fold_ix, fold in enumerate(folds):
        if symbol + "_" + str(fold_ix) not in models:
            continue
        model = models[symbol + "_" + str(fold_ix)]
        prediction = model.get_prediction(fragments, transcriptome, cell_ixs=fold["cells_test"])
        expected = prediction["expected"].values[:, 0]
        rows.append({
            "symbol": symbol,
            "fold": fold_ix,
            "cor": cor(prediction["predicted"].values[:, 0], expected),
            "cor_n_fragments": cor(prediction["n_fragments"].values[:, 0], expected),
        })
if not rows:
    sys.exit(f"no trained models found in {models_dir}")
per_fold = pd.DataFrame(rows)
per_fold["deltacor"] = per_fold["cor"] - per_fold["cor_n_fragments"]
per_fold.to_csv(out_dir / "correlations_per-fold.tsv", sep="\t", index=False)

gene_cors = per_fold.groupby("symbol", sort=False)[["cor", "cor_n_fragments", "deltacor"]].mean()
gene_cors["n_folds"] = per_fold.groupby("symbol", sort=False)["cor"].count()
gene_cors.to_csv(out_dir / "correlations.tsv", sep="\t")

# plot from https://chromatinhd.org/quickstart/3_pred/#some-quality-checks
# (genes with a NaN mean are left out; y axis not cut at 0 so negative cors show)
def correlation_plot(gene_cors, file):
    plotdata = gene_cors[["cor_n_fragments", "cor"]].dropna()
    fig, ax = plt.subplots(figsize=(4, 4))
    for name, group in plotdata.iterrows():
        ax.plot([0, 1], group[["cor_n_fragments", "cor"]], color="#3338", zorder=0, marker="o", markersize=2)
    ax.boxplot(
        plotdata.values,
        positions=[0, 1],
        widths=0.1,
        showfliers=False,
        showmeans=True,
        meanline=True,
        meanprops={"color": "red", "linewidth": 2},
    )
    ax.axhline(0, color="#999", linewidth=0.5, zorder=-1)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["# fragments", "ChromatinHD-pred"])
    ax.set_ylabel("$cor$")
    fig.savefig(file, dpi=150, bbox_inches="tight")
    plt.close(fig)

correlation_plot(gene_cors, out_dir / "correlation_plot.png")

print(f"{len(per_fold)} models ({gene_cors.shape[0]} genes), "
      f"{per_fold['cor'].isna().sum()} NaN fold cors, "
      f"{(gene_cors['n_folds'] == 0).sum()} genes without any valid fold -> {out_dir}")
