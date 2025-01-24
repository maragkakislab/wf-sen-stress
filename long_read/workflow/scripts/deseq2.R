#!/usr/bin/env Rscript

##############################################
# ORIGINALY FROM leetaiyi wf-mouse-aging
# modified by MJP

library(stringr)
library("DESeq2")
library(gtools)

####DEPRACATED: Using direct snakemake@ calls now
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


#### load in data files; explictily naming working directory may be redunandant
wd <- paste0(getwd(),"/")
print(wd)
print(snakemake@input[["counts"]]) # for debugging
print(snakemake@input[["metafile"]]) # for debugging
data_file <- read.delim(paste0(wd,snakemake@input[["counts"]]), header = TRUE, sep = snakemake@params[["delim"]], stringsAsFactors = F)
metafile <- read.delim(paste0(wd,snakemake@input[["metafile"]]), header = TRUE, sep = snakemake@params[["delim"]], stringsAsFactors = F)
print(colnames(metafile)) #for debugging
print(colnames(data_file)) # for debugging

#### prepare data_file for DESeq2
data_matrix <- as.matrix(data_file[,-1])
rownames(data_matrix)<- rownames(data_file$gene)

#### DESeq 
print(snakemake@params[["model"]])
dds <- DESeqDataSetFromMatrix(data_matrix, colData = metafile, design = as.formula(snakemake@params[["model"]]))
dds <- DESeq(dds)

####DEPRACATED use permutations now
#### iterates all comparisons from available conditions
# combinations <- combn(metafile$condition, 2, simplify = F)%>%
#   unique()

# ### function for selecting comparisons
# remove_same_same_contrasts<-function(combinations){
#   for (i in seq_along(combinations)) {
#     if (combinations[[i]][[1]] == combinations[[i]][[2]]) {
#       combinations[[i]] <- NA
#     }
#   }
#   ind<-is.na(combinations)
#   filtered<-combinations[!ind]
#   return(filtered)
# }

#### iterates all permutations of available contrasts
conditions <- unique(metafile$condition)


permutations <- permutations(n = length(conditions), 
                                  r = 2, 
                                  v = conditions)
permutations_list <- split(permutations, row(permutations))

### generates log2FC comparison for all available comparisons defined above.
print(permutations_list[[1]][[1]])

for (i in seq_along(permutations_list)) {
  df <- as.data.frame(results(dds,contrast=c("condition",permutations_list[[i]][[1]],permutations_list[[i]][[2]])))
  name <- paste(permutations_list[[i]][[1]],permutations_list[[i]][[2]],sep="_")
  file_path <- file.path(snakemake@params[["odir"]], paste0(name, ".txt"))
  write.table(df, file = file_path, row.names = TRUE, sep = snakemake@params[["delim"]])
}

### DEPRACATED WITH DIRECTORY CALL
### since output txt files will be variable based on available contrasts
### this file is used to let snakemake know the rule is completed
# comparision < - resultsNames(dds)
# write.table(comparison,file=paste0(odir,"contrasts.txt"), sep = snakemake@params[["delim"]])






