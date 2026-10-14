set -euo pipefail
RDS=/home/dsalinas/libault-group-latest-rds/integrated_multiome_annotated_minatacFALSE_15unknown_res04.rds
H5MU=/home/dsalinas/libault_multiome/integrated_multiome_annotated_minatacFALSE_15unknown_res04.h5mu
# Expression written to the RNA modality:
#   lognorm = log-normalise the raw RNA counts
#   sct     = use the object's SCTransform-normalised data
METHOD=lognorm
Rscript step-0_convert-rds-to-h5mu.R $RDS $H5MU $METHOD
