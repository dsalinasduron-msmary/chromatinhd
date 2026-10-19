# Step 2: ChromatinHD data for steps 3-5 (regions, fragments, transcriptome, folds,
# genes.tsv) in DATA (can use a lot of space).
set -euo pipefail
cd "$(dirname "$0")" # run from the repo dir (scripts, chd_patches.py, paths.sh)
source ./paths.sh # RDS, FRAGMENTS, METHOD, WINDOW and the input/output paths

python step-2_training-setup.py $H5MU $TRANSCRIPTS $FRAGMENTS $DATA $WINDOW
