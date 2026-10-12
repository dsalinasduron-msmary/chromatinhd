# Run correlation analysis. Note: this is NOT predictivity analysis.
# Correlation analysis reflects how accurately the model can predict the gene expression with __all__ 
# data available. Predictivity analysis examines the effect of hiding some data from  the model.
BATCH_SIZE=1
START_IDX=2 # start reading genes from line 2 onwards (skip header)
BATCH_FILE="batch.tsv"
GENE_FILE="test-genes.tsv"
CORRELATION_FILE_BATCH="model-window-10kb_correlations_batch.tsv"
CORRELATION_FILE="window-10kb_correlations.tsv"
TOTAL_LINES=$( wc -l $GENE_FILE | cut -d" " -f1 )

echo "symbol	cor	cor_n_fragments	deltacor	symbol" > $CORRELATION_FILE
while [ $START_IDX -le $TOTAL_LINES ]
do
echo $START_IDX of $TOTAL_LINES genes
echo GENE > $BATCH_FILE
tail -n +$START_IDX $GENE_FILE | head -n $BATCH_SIZE >> $BATCH_FILE
cat $BATCH_FILE
python3 step-4_correlation-analysis.py
if [ -f "$CORRELATION_FILE_BATCH" ]; then
	tail -n +2 $CORRELATION_FILE_BATCH >> $CORRELATION_FILE
fi
START_IDX=$(($START_IDX + $BATCH_SIZE))
done
