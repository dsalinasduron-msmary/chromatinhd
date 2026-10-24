# ChromatinHD soybean multiome pipeline

Trains [ChromatinHD](https://chromatinhd.org) models that predict each gene's expression
from the chromatin accessibility (ATAC fragments) around its transcription start site, in
soybean nodule single-cell multiome data. It then asks how well each gene can be predicted
(step 4) and which parts of the region the prediction relies on (step 5).

Input: an integrated Seurat/Signac object (`.rds`, assays RNA / SCT / ATAC with the
embedded Phytozome v13 gene annotation) and its ATAC fragments file.

## Pipeline

Each step has a shell wrapper `step-N.sh` (no arguments) that calls an R or Python
script with command-line arguments. Run them in order:

| Step | Wrapper → script | Reads | Writes |
|---|---|---|---|
| 0 | `step-0.sh` → `step-0_convert-rds-to-h5mu.R` (R) | `RDS` | `H5MU`: RNA expression (per `METHOD`) + ATAC peaks |
| 1a | `step-1.sh` → `step-1_annotation-to-transcripts.R` (R) | `RDS` | `CANDIDATES`: every annotated transcript of every RNA gene |
| 1b | `step-1.sh` → `step-1_select-tss.py` | `CANDIDATES`, `FRAGMENTS` | `TRANSCRIPTS`: one TSS per gene (most fragments within ±100 bp; ties → longest transcript) |
| 2 | `step-2.sh` → `step-2_training-setup.py` | `H5MU`, `TRANSCRIPTS`, `FRAGMENTS` | `DATA/`: regions (TSS ± `WINDOW`), fragments, transcriptome, 5 cell folds, `genes.tsv` |
| 3 | `step-3.sh` → `step-3_train-models.py` | `DATA/`, gene list | `DATA/models/`: one model per gene × fold (GPU) |
| 4 | `step-4.sh` → `step-4_correlation-analysis.py` | `DATA/`, models | `DATA/correlations/` |
| 5 | `step-5.sh` → `step-5_predictivity-analysis.py` | `DATA/`, models, step-4 table | `DATA/predictivity/` |

Usage of each script is in its header comment (`Usage: ...`).

## Configuration: `paths.sh`

All wrappers source `paths.sh`; edit it, not the wrappers, for anything shared:

- `RDS`, `FRAGMENTS`: inputs. The fragments file must be bgzipped and tabix-indexed
  (`atac_fragments.tsv.gz.tbi` next to it; make one with `tabix -p bed <file>`).
- `METHOD`: expression to model. `lognorm` = log-normalised raw RNA counts;
  `sct` = the object's SCTransform data.
- `WINDOW`: region per gene = TSS ± `WINDOW` bp (default 10000). ChromatinHD's default
  fragment encoder bins (100–5000 bp) only fit multiples of 10 kb; for other windows step 3
  prints a warning and uses custom bin widths that divide the window (e.g. 3000 →
  100, 200, 500, 1000, 1500, 3000 bp).
- `H5MU_DIR`, `WORK`: where outputs go.

Output names are derived from these (`<rds name>_<METHOD>.h5mu`,
`WORK/<METHOD>_window-<WINDOW>bp/`), so runs with a different `METHOD` or `WINDOW` sit
side by side instead of overwriting each other. `TRANSCRIPTS` (step 1) depends on neither.

Settings for a single step stay in its wrapper:

- `step-3.sh`: `GENE_FILE` (genes to train; default `DATA/genes.tsv` = all genes from
  step 2; point it at a shorter list, column `GENE`, to train fewer), `BATCH_SIZE`
  (genes per Python run; batching avoids "too many open files").
- `step-5.sh`: `MIN_COR`, `MIN_DELTACOR` (which genes to analyse, see step 5 below).

The wrappers can be started from any directory (they `cd` to the repo first).

## Outputs

**Step 3** — `DATA/models/`. Reruns skip models that already exist, so an interrupted
run can be restarted, and genes can be added later.

**Step 4** — `DATA/correlations/`, scored on each fold's held-out test cells:

- `correlations_per-fold.tsv`: `symbol, fold, cor, cor_n_fragments, deltacor` per gene ×
  fold. `cor` = correlation of the ChromatinHD prediction with expression;
  `cor_n_fragments` = correlation of the fragment count in the region with expression
  (baseline); `deltacor = cor - cor_n_fragments`. NaN where expression (or the
  prediction) is constant across a fold's test cells.
- `correlations.tsv`: one row per gene, means over the folds that are not NaN, and
  `n_folds` (how many folds counted).
- `correlation_plot.png`: per-gene `cor_n_fragments` vs `cor`.

**Step 5** — `DATA/predictivity/`, for genes with `cor >= MIN_COR` and
`deltacor >= MIN_DELTACOR` in `correlations.tsv` (default 0.2 / 0.05). Every model hides
windows of 50–500 bp one at a time and records the drop in prediction quality.

- `selected_genes.tsv`: the genes analysed.
- `predictivity.tsv`: per gene every 10 bp: `position` (bp from the TSS; + = downstream
  on the gene's strand), `deltacor` (< 0: hiding this window makes the prediction worse),
  `lost` (fragments hidden per cell), `effect` (> 0: accessibility here goes with higher
  expression).
- `predictive_windows.tsv`: stretches with `deltacor < -0.001` (ChromatinHD default),
  by `effect_direction` (+1 / -1), with TSS-relative `start`/`end` and genome
  coordinates `chrom`, `start_genome`, `end_genome`.
- `plots/<gene>.png`: fragment pileup and predictivity along the region.
- `regionmultiwindow/`: ChromatinHD's stored results; reruns skip genes already scored.

## Environment

- **R 4.5** with Seurat, Signac, MuDataSeurat, Rsamtools, ... (steps 0, 1a). Versions
  and install notes: `packages_R.txt`.
- **Python 3.12** with chromatinhd 0.4.3, torch 2.14.1 (CUDA 13 build) (steps 1b–5).
  Exact pins: `requirements.txt` (compiled from `requirements.in`; recompile command in
  its header). Notes: `packages_Python.txt`.
- `chd-soybean.def` (Singularity/Apptainer) builds both; `Dockerfile` has the same
  steps for Docker (keep the two in sync). The Singularity image
  was built and tested on the production server (2026-10-20); the Dockerfile has not been
  built.

### Running on the production server (Singularity/Apptainer)

Build the image once (needs internet access; ~30–60 min, mostly compiling R packages):

```bash
bash build-image.sh
```

It writes `$WORK/chd-soybean.sif` (`WORK` from `paths.sh`) and keeps its temporary
build files (~15–20 GB at the peak) and download cache under `WORK` instead of `/tmp`
and your home directory. It builds with `--fakeroot`; if the admins have not enabled
that, set `BUILD_OPTS=""` in the script and run it with `sudo`. The build ends with a
test that loads the R and Python packages.

Then run each step from the repo directory (the scripts are not in the image):

```bash
singularity exec --nv --bind /falafel /falafel/cornbert/chd-soybean.sif bash step-0.sh
singularity exec --nv --bind /falafel /falafel/cornbert/chd-soybean.sif bash step-1.sh
# ... through step-5.sh
```

- `--nv` gives the container the GPU (steps 3–5).
- `--bind` makes directories outside your home visible inside the container: add every
  directory `paths.sh` refers to (e.g. `--bind /falafel,/cluster`) unless the site
  already mounts them.
- GPU: the CUDA 13 build of torch supports compute capability 7.5 and newer (T4,
  RTX 20xx and later, A100, H100, ...), not P100 / V100, and needs an NVIDIA driver
  ≥ 580 on the host (check with `nvidia-smi`). Without a GPU, ChromatinHD falls back to
  the CPU (not tested; expect it to be much slower).

### chromatinhd 0.4.3 workarounds

- `chd_patches.py`, imported by steps 3–5 right after chromatinhd, fixes two upstream
  bugs: the pred models call `torch_scatter` although its import was removed (replaced
  with pure-torch equivalents), and the step-5 plot panels (polyptich `Panel.ax` never
  initialised). Remove when chromatinhd fixes them.
- Pinned `anndata<0.13` (0.13 breaks `Transcriptome.from_adata`) and `zarr<3`
  (chromatinhd's sparse result storage, used in step 5, is written for zarr 2).

## Resources

Measured on an 870-cell, 500-gene test subset (RTX 2070 Super); the full data
(11,827 cells, ~40,000 genes) will be larger and slower:

- Step 3: ~10 s per gene (5 models); per model ~0.44 MB on disk and 8 files.
  All genes ≈ 200,000 models ≈ 100 GB and 1.6 million files (inodes) — check the
  server's disk and file-count quotas first, or train a gene subset.
- Step 4: ~47 s for 2,500 models.
- Step 5: ~1.5 min for 73 genes.

## Testing

`make-test-subset.R` builds a small dataset (a few clusters, the most variable genes,
peaks and fragments near them) from the full object:

```bash
Rscript make-test-subset.R <full.rds> <fragments.tsv.gz> test-data/test-subset_clusters-10-15-16_genes-500 \
    500 50000 10 15 16
```

Point `paths.sh` at the resulting `.rds` / `_fragments.tsv.gz` (and `H5MU_DIR` / `WORK`
at a scratch directory) to run the whole pipeline on it in minutes; shorten
`GENE_FILE` in `step-3.sh` to train only a few genes.

`inspect.R <object.rds>` prints an object's assays, layers, metadata, clusters,
barcodes and fragment paths — useful to check an object (e.g. on the production
server) before step 0.
