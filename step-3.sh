BATCH_SIZE=500
START_IDX=2 # start reading genes from line 2 onwards (skip header)
BATCH_FILE="genes-of-interest_batch.tsv"
GENE_FILE="test-genes.tsv"
TOTAL_LINES=$( wc -l $GENE_FILE | cut -d" " -f1 )

while [ $START_IDX -le $TOTAL_LINES ]
do
echo $START_IDX of $TOTAL_LINES genes
echo GENE > $BATCH_FILE
tail -n +$START_IDX $GENE_FILE | head -n $BATCH_SIZE >> $BATCH_FILE
python3 step-3_train-models.py
START_IDX=$(($START_IDX + $BATCH_SIZE))
done
