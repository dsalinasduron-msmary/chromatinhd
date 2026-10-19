# Step 0: Seurat/Signac .rds -> .h5mu (RNA expression per METHOD, ATAC peaks).
set -euo pipefail
cd "$(dirname "$0")" # run from the repo dir (scripts, chd_patches.py, paths.sh)
source ./paths.sh # RDS, FRAGMENTS, METHOD, WINDOW and the input/output paths

mkdir -p $H5MU_DIR
Rscript step-0_convert-rds-to-h5mu.R $RDS $H5MU $METHOD
