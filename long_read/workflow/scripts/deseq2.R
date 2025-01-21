#!/usr/bin/env Rscript

##############################################
# ORIGINALY FROM leetaiyi wf-mouse-aging
# modified by MJP


suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
  library(tidyverse)
})

option_list <- list(
  make_option(c("-i","--data"),
              help="Path to aggregate datafile used in comparsion;\ncolumns:gene, sample_n,[addtional samples], ..."),
  make_option(c("-m","--metadata"),
              help="Path to desq2 formated metadata file;\ncolumns: sample, condition, batch, [additional columns]"),
  make_option(c("-r","--model", default = "~ condition"),
              help="Model to use for DEseq2 [default: ~ condition]"),
  make_option(c("-d","--delim", default = "\t"),
              help="Delimiter [default tab]"),
  make_option(c("-o","--odir"), default = ".",
              help="Path to output directory [default current directory]")
)

opt <- parse_args(OptionParser(option_list = option_list))

data_file <- read.delim(opt$data, header = TRUE, sep = opt$delim, stringsAsFactors = F)
metadata <- read.delim(opt$metadata, header = TRUE, sep = opt$delim, stringsAsFactors = T)

#Check whether metadata rownames exist in count data
if (!identical(colnames(metadata[,-1]),data_file$sample)) {
  stop("metadata samples do not match order/content of data_file sample colummns")
}

#Prepare data_file for DESeq2
data_matrix <- as.matrix(data_file[,-1])
rownames(data_matrix)<- data_file$gene


#DESeq 

dds <- DESeqDataSetFromMatrix(data_matrix, colData = metadata, design = opt$model)
dds <- DESeq(dds)

#Generate contrasts
odir <- opt$odir

for (i in resultsNames(dds)) {
  df <- as.data.frame(results(dds, name = i))
  file_path <- file.path(output_dir, paste0(i, ".txt"))
  write.table(df, file = file_path, row.names = TRUE, sep = opt$delim)
}

### Since output txt files will be variable based on available contrasts
### This file is used to let snakemake know the rule is completed
write.table(comparison,file=paste0(odir,"contrasts.txt"), sep = opt$delim)






