# Run predictivity analysis: for each selected gene, which windows of the region the model relies on.
# Genes are selected from the step-4 per-gene correlations by the two thresholds below.
set -euo pipefail

DATA=/falafel/cornbert/window-10kb # step-2 output (same as DATA in step-3.sh / step-4.sh)
MODELS=$DATA/models # step-3 output (same as MODELS in step-3.sh / step-4.sh)
CORRELATIONS=$DATA/correlations/correlations.tsv # step-4 per-gene table
MIN_COR=0.2 # keep genes whose ChromatinHD prediction has cor >= MIN_COR
MIN_DELTACOR=0.05 # ... and beats counting fragments by >= MIN_DELTACOR
OUT=$DATA/predictivity # gets selected_genes.tsv, predictivity.tsv, predictive_windows.tsv, plots/

python step-5_predictivity-analysis.py $DATA $MODELS $CORRELATIONS $MIN_COR $MIN_DELTACOR $OUT
