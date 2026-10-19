# Run correlation analysis. Note: this is NOT predictivity analysis.
# Correlation analysis reflects how accurately the model can predict the gene expression with __all__
# data available. Predictivity analysis (step 5) examines the effect of hiding some data from the model.
# Scores every model trained in MODELS (one pass; fast enough that no batching is needed).
set -euo pipefail
cd "$(dirname "$0")" # run from the repo dir (scripts, chd_patches.py, paths.sh)
source ./paths.sh # RDS, FRAGMENTS, METHOD, WINDOW and the input/output paths

OUT=$DATA/correlations # gets correlations_per-fold.tsv, correlations.tsv, correlation_plot.png

python step-4_correlation-analysis.py $DATA $MODELS $OUT
