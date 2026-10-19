# Run predictivity analysis: for each selected gene, which windows of the region the model relies on.
# Genes are selected from the step-4 per-gene correlations by the two thresholds below.
set -euo pipefail
cd "$(dirname "$0")" # run from the repo dir (scripts, chd_patches.py, paths.sh)
source ./paths.sh # RDS, FRAGMENTS, METHOD, WINDOW and the input/output paths

CORRELATIONS=$DATA/correlations/correlations.tsv # step-4 per-gene table
MIN_COR=0.2 # keep genes whose ChromatinHD prediction has cor >= MIN_COR
MIN_DELTACOR=0.05 # ... and beats counting fragments by >= MIN_DELTACOR
OUT=$DATA/predictivity # gets selected_genes.tsv, predictivity.tsv, predictive_windows.tsv, plots/

python step-5_predictivity-analysis.py $DATA $MODELS $CORRELATIONS $MIN_COR $MIN_DELTACOR $OUT
