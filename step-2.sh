set -euo pipefail
# step 0 output (expression to model in mod/RNA .X)
H5MU=/home/dsalinas/libault_multiome/integrated_multiome_annotated_minatacFALSE_15unknown_res04.h5mu
# step 1 output
TRANSCRIPTS=/falafel/cornbert/transcripts.csv
# bgzipped and tabix-indexed (.tbi)
FRAGMENTS=/falafel/pgrp_singlecell/soybean_nodule_multiome-redo/atac_fragments.tsv.gz
# ChromatinHD data for steps 3-5 (can use a lot of space); name it after the setup
OUT=/falafel/cornbert/window-10kb
WINDOW=10000
python step-2_training-setup.py $H5MU $TRANSCRIPTS $FRAGMENTS $OUT $WINDOW
