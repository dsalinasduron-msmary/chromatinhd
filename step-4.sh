# Run correlation analysis. Note: this is NOT predictivity analysis.
# Correlation analysis reflects how accurately the model can predict the gene expression with __all__
# data available. Predictivity analysis (step 5) examines the effect of hiding some data from the model.
# Scores every model trained in MODELS (one pass; fast enough that no batching is needed).
set -euo pipefail

DATA=/falafel/cornbert/window-10kb # step-2 output (same as DATA in step-3.sh)
MODELS=$DATA/models # step-3 output (same as MODELS in step-3.sh)
OUT=$DATA/correlations # gets correlations_per-fold.tsv, correlations.tsv, correlation_plot.png

python step-4_correlation-analysis.py $DATA $MODELS $OUT
