# Step 1: one TSS per gene from the object's gene annotation.
#   1a: all candidate transcripts (-> CANDIDATES)
#   1b: per gene the transcript whose TSS has most fragments (-> TRANSCRIPTS)
set -euo pipefail
cd "$(dirname "$0")" # run from the repo dir (scripts, chd_patches.py, paths.sh)
source ./paths.sh # RDS, FRAGMENTS, METHOD, WINDOW and the input/output paths

mkdir -p $WORK
Rscript step-1_annotation-to-transcripts.R $RDS $CANDIDATES
python step-1_select-tss.py $CANDIDATES $FRAGMENTS $TRANSCRIPTS
