set -euo pipefail
RDS=/home/dsalinas/libault-group-latest-rds/integrated_multiome_annotated_minatacFALSE_15unknown_res04.rds
# bgzipped and tabix-indexed (.tbi)
FRAGMENTS=/falafel/pgrp_singlecell/soybean_nodule_multiome-redo/atac_fragments.tsv.gz
# all candidate transcripts (step 1a) and one TSS per gene for step 2 (step 1b)
CANDIDATES=/falafel/cornbert/transcript-candidates.csv
TRANSCRIPTS=/falafel/cornbert/transcripts.csv
Rscript step-1_annotation-to-transcripts.R $RDS $CANDIDATES
python step-1_select-tss.py $CANDIDATES $FRAGMENTS $TRANSCRIPTS
