library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(readr)

pt <- "/falafel/cornbert"
fn <- paste(pt,"glyma.Wm82.gnm4.ann1.T8TQ.gene_models_main.bed",sep="/")
df <- read_tsv(fn,col_names=F)

# X1: chr "glyma.Wm82.gnm4.Gm01" --> chrom
# X2: num 27655 --> start, ---> tss
# X3: num 28218 --> end
# X4: chr "glyma.Wm82.gnm4.ann1.Glyma.01G000100.1" --> ensembl_gene_id.v
# X5: num 0  [DROP]
# X6: chr "-" --> strand
# X7: chr "glyma.Wm82.gnm4.ann1.Glyma.01G000100" --> ensembl_gene_id

#g1 <- "glyma.Wm82.gnm4.ann1.Glyma.19G087000" has multiple versions
#g2 <- "glyma.Wm82.gnm4.ann1.Glyma.20G229000" has multiple versions
transcript.df <- df %>%
	           rename(chrom=X1) %>%
	             rename(start=X2) %>%
	               mutate(tss=start) %>%
	                 rename(end=X3) %>%
	                   rename(strand=X6) %>%
	                     rename(ensembl_gene_id=X7) %>%
	                       rename(transcript=X4) %>% 
	                         select(-X5) %>% 
				   group_by(ensembl_gene_id,chrom,strand) %>%
				       slice(1) %>%
				         ungroup()
## needed 
# ensembl_gene_id start     end         strand symbol   transcript_biotype chrom  tss
# GLYMA_01G143900 47782788  47787690       1   LBD3     protein_coding     chr1   47782788
transcript.df

## current
#ensembl_gene_id                      chrom        strand  start    end    tss
#glyma.Wm82.gnm4.ann1.Glyma.01G000100 glyma.Wm82.~ -       27655  28218  27655
#glyma.Wm82.gnm4.ann1.Glyma.01G000137 glyma.Wm82.~ -       52484  52601  52484

f.1 <- function(x) {
	prefix <-  "glyma.Wm82.gnm4.ann1."
	return (x %>% str_replace(prefix,"") %>% str_replace("[.]","_") %>% str_to_upper())
}
f.2 <- function(x) {
	prefix <-  "glyma.Wm82.gnm4."
	return (x %>% str_replace(prefix,"") )
}
f.3 <- function(x) {
	prefix <- "glyma.Wm82.gnm4.ann1."
	return ( x %>% str_replace(prefix,"") %>% str_replace("[.]","_") %>% str_to_upper() )
}
transcript.df %<>%
	mutate(ensembl_gene_id = f.1(ensembl_gene_id) ) %>%
	  mutate(chrom = f.2(chrom) ) %>%
	    mutate(strand = ifelse(strand == "+",1,-1) ) %>%
	      mutate(symbol=ensembl_gene_id) %>%
	        mutate(transcript_biotype="protein_coding") %>%
		  mutate(transcript = f.3(transcript) )

out.file <- paste(pt,"transcripts.csv",sep="/")
transcript.df %>% write_csv(out.file)
