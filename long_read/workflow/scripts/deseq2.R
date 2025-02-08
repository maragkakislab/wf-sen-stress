#!/usr/bin/env Rscript

##############################################
# ORIGINALY FROM leetaiyi wf-mouse-aging
# modified by MJP

library(stringr)
library("DESeq2")
library(gtools)

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
rownames(data_matrix)<- data_file$gene

#### DESeq 
print(snakemake@params[["model"]])
dds <- DESeqDataSetFromMatrix(data_matrix, colData = metafile, design = as.formula(snakemake@params[["model"]]))
dds <- DESeq(dds)

#### iterates all permutations of available contrasts
conditions <- unique(metafile$condition)


permutations <- permutations(n = length(conditions), 
                                  r = 2, 
                                  v = conditions)
permutations_list <- split(permutations, row(permutations))

### generates log2FC comparison for all available comparisons defined above.
print(permutations_list[[1]][[1]])
id<-read.delim(snakemake@input[["ENSG_metadata"]], header = TRUE, sep = snakemake@params[["delim"]], stringsAsFactors = F)
if (!dir.exists(snakemake@params[["odir"]])) {
  dir.create(snakemake@params[["odir"]], recursive = TRUE)
}

for (i in seq_along(permutations_list)) {
  df<- as.data.frame(results(dds,contrast=c("condition",permutations_list[[i]][[1]],permutations_list[[i]][[2]])))
  ensembl<-rownames(df)
  df$ensembl<-ensembl
  annotated_df<-merge(id, df, by.y= "ensembl", by.x="ensembl_gene_id")
  name <- paste(permutations_list[[i]][[1]],permutations_list[[i]][[2]],sep="_")
  file_path <- file.path(snakemake@params[["odir"]], paste0(name, ".dge.txt"))
  print(file_path)
  write.table(annotated_df, file = file_path, row.names = FALSE, sep = snakemake@params[["delim"]])
}

### since output txt files will be variable based on available contrasts
### this file is used to let snakemake know the rule is completed
write.table(permutations,file=snakemake@output[["permutations"]], row.names=FALSE, sep = snakemake@params[["delim"]])






