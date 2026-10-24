# Shared settings for step-0.sh .. step-5.sh (each wrapper sources this file).
# Edit inputs and choices here; step-specific settings (BATCH_SIZE, thresholds, ...)
# stay in the individual wrappers.

# --- inputs ---
# Seurat/Signac object (RNA, SCT, ATAC with the gene annotation)
RDS=/home/dsalinas/libault-group-latest-rds/integrated_multiome_annotated_links_doublet2p4_cluster_reordered.rds

# ATAC fragments, bgzipped and tabix-indexed (.tbi next to it)
FRAGMENTS=/falafel/cornbert/atac_fragments.tsv.gz

# --- choices ---
# expression written to the RNA modality in step 0:
#   lognorm = log-normalise the raw RNA counts
#   sct     = use the object's SCTransform-normalised data
METHOD=sct
# region modelled per gene: TSS +/- WINDOW bp (step 2). Not a multiple of 10000 ->
# step 3 warns and uses custom encoder bin widths that fit the window
WINDOW=3000

# --- output locations ---
H5MU_DIR=/falafel/cornbert/chromatinhd_artifacts # step 0 output
WORK=/falafel/cornbert/chromatinhd_claude_artifacts/work # steps 1-5 output (can use a lot of space)

# --- derived names (no need to edit): METHOD and WINDOW are part of the names so
# runs with other settings do not overwrite each other ---
H5MU=$H5MU_DIR/$(basename $RDS .rds)_$METHOD.h5mu # step 0 -> step 2
CANDIDATES=$WORK/transcript-candidates.csv # step 1a -> step 1b
TRANSCRIPTS=$WORK/transcripts.csv # step 1 -> step 2 (does not depend on METHOD/WINDOW)
DATA=$WORK/${METHOD}_window-${WINDOW}bp # step 2 -> steps 3-5
MODELS=$DATA/models # step 3 -> steps 4-5
