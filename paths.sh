# Shared settings for step-0.sh .. step-5.sh (each wrapper sources this file).
# Edit inputs and choices here; step-specific settings (BATCH_SIZE, thresholds, ...)
# stay in the individual wrappers.

# --- inputs ---
# Seurat/Signac object (RNA, SCT, ATAC with the gene annotation)
RDS=/home/dsalinas/libault-group-latest-rds/integrated_multiome_annotated_minatacFALSE_15unknown_res04.rds
# ATAC fragments, bgzipped and tabix-indexed (.tbi next to it)
FRAGMENTS=/falafel/pgrp_singlecell/soybean_nodule_multiome-redo/atac_fragments.tsv.gz

# --- choices ---
# expression written to the RNA modality in step 0:
#   lognorm = log-normalise the raw RNA counts
#   sct     = use the object's SCTransform-normalised data
METHOD=lognorm
# region modelled per gene: TSS +/- WINDOW bp (step 2)
WINDOW=10000

# --- output locations ---
H5MU_DIR=/home/dsalinas/libault_multiome # step 0 output
WORK=/falafel/cornbert # steps 1-5 output (can use a lot of space)

# --- derived names (no need to edit): METHOD and WINDOW are part of the names so
# runs with other settings do not overwrite each other ---
H5MU=$H5MU_DIR/$(basename $RDS .rds)_$METHOD.h5mu # step 0 -> step 2
CANDIDATES=$WORK/transcript-candidates.csv # step 1a -> step 1b
TRANSCRIPTS=$WORK/transcripts.csv # step 1 -> step 2 (does not depend on METHOD/WINDOW)
DATA=$WORK/${METHOD}_window-${WINDOW}bp # step 2 -> steps 3-5
MODELS=$DATA/models # step 3 -> steps 4-5
