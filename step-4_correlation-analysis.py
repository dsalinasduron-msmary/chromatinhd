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
model_label = "homoelogs" # added as suffix to the training artifacts
model_list = "batch.tsv"

# recover models
# do this first so we can examine only regions that have models
models = chd.models.pred.model.multiscale.Models(pt + "/" + "models_" + model_label, reset=False)

# recover regions and select custom subset of models to be examined (regions_oi )
regions = chd.data.regions.Regions(pt + "/" + "regions_" + model_label)
regions_oi = list(pd.read_csv(model_list)["GENE"])

# recover transcriptome
transcriptome = chd.data.Transcriptome(pt + "/" + "transcriptome_" + model_label)

# recover folds
folds = chd.data.folds.Folds(pt + "/" + "folds_" + model_label + "/5x1") 

# recover fragments
fragments = chd.data.Fragments(pt + "/" + "fragments_" + model_label)
fragments.create_regionxcell_indptr()

# filter to focus on subset of genes (models of regions_oi)
f = fragments.filter_regions(regions.filter(regions_oi))
f.create_regionxcell_indptr()
f.obs["cells"] = f.obs.index
f.obs.set_index("cells",inplace=True)
t = transcriptome.filter_genes(regions_oi)

gene_cors = models.get_region_cors(f, t, folds)
gene_cors.dropna(inplace=True)
#gene_cors["symbol"] = gene_cors.index
#gene_cors.to_csv(model_label + "_correlations_batch.tsv",sep="\t")

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

