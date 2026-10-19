set -euo pipefail
# ChromatinHD data from step 2
DATA=/falafel/cornbert/window-10kb
# genes to train, one per line under a GENE header (GLYMA_XXGXXXXXX names)
GENE_FILE=test-genes.tsv
# trained models (can use a lot of space)
MODELS=$DATA/models
# genes per run of the Python script (batches avoid a "too many open files" error)
BATCH_SIZE=500
BATCH_FILE=$(mktemp --suffix=.tsv)
trap 'rm -f $BATCH_FILE' EXIT

START_IDX=2 # start reading genes from line 2 onwards (skip header)
TOTAL_LINES=$( wc -l < $GENE_FILE )
while [ $START_IDX -le $TOTAL_LINES ]
do
END_IDX=$(( $START_IDX + $BATCH_SIZE - 1 < $TOTAL_LINES ? $START_IDX + $BATCH_SIZE - 1 : $TOTAL_LINES ))
echo genes $(($START_IDX - 1))-$(($END_IDX - 1)) of $(($TOTAL_LINES - 1))
echo GENE > $BATCH_FILE
sed -n "${START_IDX},${END_IDX}p" $GENE_FILE >> $BATCH_FILE
python step-3_train-models.py $DATA $BATCH_FILE $MODELS
START_IDX=$(($START_IDX + $BATCH_SIZE))
done
