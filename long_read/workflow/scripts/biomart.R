library(biomaRt)
library(tidyverse)
ensembl_to_gene<-function(df,mart){
  id<-getBM(attributes = c(snakemake@params[["identifier"]],'description','gene_biotype','ensembl_gene_id'),
          filters = 'ensembl_gene_id',
          values = ensembl_gene, 
          mart = mart)
  return(id)
}

mart <- useEnsembl(biomart = "genes", dataset = snakemake@params[["genome"]])
ensembl_table <- read.delim(snakemake@input[["transcript_tab"]],sep="\t",header=TRUE)
ensembl_gene<-ensembl_table$gene %>% unique()

write.table(ensembl_to_gene(ensembl_gene,mart), file = snakemake@output[["ENSG_metadata"]], row.names = FALSE, sep = "\t")