#!/usr/bin/env Rscript

##############################################
# ORIGINALY FROM leetaiyi wf-mouse-aging
# modified by MJP

library(stringr)
library("DESeq2")

# option_list <- list(
#   make_option(c("-i","--data"),
#               help="Path to aggregate datafile used in comparsion;\ncolumns:gene, sample_n,[addtional samples], ..."),
#   make_option(c("-m","--metadata"),
#               help="Path to desq2 formated metadata file;\ncolumns: sample, condition, batch, [additional columns]"),
#   make_option(c("-r","--model", default = "~ condition"),
#               help="Model to use for DEseq2 [default: ~ condition]"),
#   make_option(c("-d","--delim", default = "\t"),
#               help="Delimiter [default tab]"),
#   make_option(c("-o","--odir"), default = ".",
#               help="Path to output directory [default current directory]")
# )

# opt <- parse_args(OptionParser(option_list = option_list))
wd <- paste0(getwd(),"/")
print(wd)
print(snakemake@input[["counts"]]) # for debugging
print(snakemake@input[["metafile"]]) # for debugging
data_file <- read.delim(paste0(wd,snakemake@input[["counts"]]), header = TRUE, sep = snakemake@params[["delim"]], stringsAsFactors = F)
metafile <- read.delim(paste0(wd,snakemake@input[["metafile"]]), header = TRUE, sep = snakemake@params[["delim"]], stringsAsFactors = T)

#Check whether metadata rownames exist in count data
print(colnames(metafile)) #for debugging
print(colnames(data_file)) # for debugging


# if (!identical(colnames(data_file[,-1]),metafile$sample)) {
#   stop("metafile samples do not match order/content of data_file sample colummns")
# }


#Prepare data_file for DESeq2
data_matrix <- as.matrix(data_file[,-1])
rownames(data_matrix)<- rownames(data_file$gene)

#DESeq 

dds <- DESeqDataSetFromMatrix(data_matrix, colData = metafile, design = snakemake@params[["model"]])
dds <- DESeq(dds)

#Generate contrasts

for (i in resultsNames(dds)) {
  df <- as.data.frame(results(dds, name = i))
  file_path <- file.path(wd, paste0(i, ".txt"))
  write.table(df, file = file_path, row.names = TRUE, sep = snakemake@params[["delim"]])
}

### Since output txt files will be variable based on available contrasts
### This file is used to let snakemake know the rule is completed
write.table(comparison,file=paste0(odir,"contrasts.txt"), sep = snakemake@params[["delim"]])






