import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import polyptich
import torch
import anndata as ad
import chromatinhd as chd
import chd_patches  # fixes for chromatinhd 0.4.3 (torch_scatter); see chd_patches.py
import mudata
import pandas as pd

# main folder (this is decided in training)
pt = "/falafel/cornbert" # path to put training data and models
model_label = "homoelogs" # suffix of the model folder (set in step 3)
model_list = "high-corr_high-deltacor_genes.tsv"
output_prefix = "feb-9" # attached to spreadsheets
#df = pd.read_csv("pairs.tsv",sep="\t")

# recover models
# do this first so we can examine only regions that have models
models = chd.models.pred.model.multiscale.Models(pt + "/" + "models_" + model_label, reset=False)

# recover regions and select custom subset of models to be examined (regions_oi )
regions = chd.data.regions.Regions(pt + "/" + "regions")
regions_oi = list(pd.read_csv(model_list,sep="\t")["symbol"])
#regions_oi = list(df.iloc[0])

# recover transcriptome
transcriptome = chd.data.Transcriptome(pt + "/" + "transcriptome")

# recover folds
folds = chd.data.folds.Folds(pt + "/" + "folds/5x1") 

# recover fragments
fragments = chd.data.Fragments(pt + "/" + "fragments")

# filter to focus on subset of genes (models of regions_oi)
f = fragments.filter_regions(regions.filter(regions_oi))
f.create_regionxcell_indptr()
f.obs["cells"] = f.obs.index
f.obs.set_index("cells",inplace=True)
t = transcriptome.filter_genes(regions_oi)

"""
gene_cors = models.get_region_cors(f, t, folds)
gene_cors.dropna(inplace=True)
gene_cors["symbol"] = gene_cors.index
gene_cors.to_csv(model_label + "_correlations_batch.tsv",sep="\t")
"""

# plot from https://chromatinhd.org/quickstart/3_pred/#some-quality-checks
def correlation_plot():
    fig, ax = plt.subplots(figsize=(4, 4))
    for name, group in gene_cors.iterrows():
        ax.plot([0, 1], group[["cor_n_fragments", "cor"]], color="#3338", zorder=0, marker="o", markersize=2)
    ax.boxplot(
        gene_cors[["cor_n_fragments", "cor"]].values,
        positions=[0, 1],
        widths=0.1,
        showfliers=False,
        showmeans=True,
        meanline=True,
        meanprops={"color": "red", "linewidth": 2},
    )
    ax.set_xticks([0, 1])
    ax.set_ylim(0)
    ax.set_xticklabels(["# fragments", "ChromatinHD-pred"])
    ax.set_ylabel("$cor$")
    plt.savefig("correlation_plot.png")

#correlation_plot()
# plot can be messy. skip for now.

    
censorer = chd.models.pred.interpret.MultiWindowCensorer(f.regions.window)
censorer.__class__ = chd.models.pred.interpret.censorers.MultiWindowCensorer
regionmultiwindow = chd.models.pred.interpret.RegionMultiWindow.create(
    path = models.path / "interpret" / ("regionmultiwindow_" + output_prefix),
    folds = folds,
    transcriptome = t,
    censorer = censorer,
    fragments = f,
    overwrite = False,
)

regionmultiwindow.score(
    models = models,
    regions = regions_oi,
    folds = folds,
)
regionmultiwindow.interpolate()


for symbol in regions_oi :
    fig = polyptich.grid.Figure(polyptich.grid.Grid(padding_height=0.05))
    plotdata = regionmultiwindow.get_plotdata(symbol).reset_index()
    plotdata.to_csv(output_prefix + "_" + symbol + ".csv")
    width = 10
    region = fragments.regions.coordinates.loc[symbol]
    #panel_genes = chd.plot.genome.genes.Genes.from_region(region, width=width, genome = "GRCh38")
    #fig.main.add_under(panel_genes)
    panel_pileup = chd.models.pred.plot.Pileup.from_regionmultiwindow(
        regionmultiwindow, symbol, width=width
    )
    fig.main.add_under(panel_pileup)
    panel_predictivity = chd.models.pred.plot.Predictivity.from_regionmultiwindow(
        regionmultiwindow, symbol, width=width
    )
    fig.main.add_under(panel_predictivity)
    fig.plot()
    plt.savefig(model_label + "_predictivity-plot_" + symbol +  ".png")
    #plt.savefig(model_label + "_predictivity-plot_" + "test" +  ".png")
    plt.close()

"""
df = regionmultiwindow.extract_predictive_windows()
df.to_csv(model_label + "_predictivity-window.tsv",sep="\t",index=False)
"""
