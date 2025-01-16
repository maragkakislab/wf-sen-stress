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
  make_option(c("-p", "--pcvars", default = NULL),
              help="PC plot variables, separated by commas. Example: sex,time. If left empty, PCA is not performed."),
  make_option(c("-o","--odir"), default = ".",
              help="Path to output directory [default current directory]")
)

opt <- parse_args(OptionParser(option_list = option_list))

data_file <- read.delim(opt$data, header = TRUE, sep = delim, stringsAsFactors = F)
metadata <- read.delim(opt$metadata, header = TRUE, sep = delim, stringsAsFactors = T)

#Check whether metadata rownames exist in count data
if (!identical(colnames(metadata[,-1]),data_file$sample)) {
  stop("metadata samples do not match order/content of data_file sample colummns")
}

#Prepare data_file for DESeq2
data_matrix <- as.matrix(data_file[,-1])
rownames(data_matrix)<- data_file$gene


#DESeq 

dds = DESeqDataSetFromMatrix(data_matrix, colData = metadata, design = opt$model)
dds = estimateSizeFactors(dds)
counts_norm = counts(dds, normalized = T)

#Generate contrasts






if (!is.null(model)) {
  design = as.formula(model)
  
  v = head(strsplit(model, "\\~|\\+")[[1]],1) #Find first variable in model
  t = table(metadata[,v])
  if (length(t) > 2 ) {
    warning("More than two levels of last model variable!")
  }
  if (!(baseline %in% names(t))) {
    stop("Baseline not a level in last model variable!")
  }
  notbaseline = names(t)[-which(names(t)==baseline)]
  
  dds <- DESeqDataSetFromMatrix(countData = counts_wide,
                                colData = metadata,
                                design = design)
  
  dds <- DESeq(dds)
  outres <- results(dds, contrast = c(v,notbaseline,baseline))
  outres <- data.frame(rownames(outres),outres)
  colnames(outres)[1] = ctg_col
  outres$geneid = gprofiler2::gconvert(rownames(outres),
                                       organism = species, target = "ENSG",filter_na = F)$name
  
  g = ggplot(outres, aes(x = log2FoldChange, y = -log10(padj))) + geom_point() +
    geom_hline(yintercept=-log10(0.05), color = "red", size = 1) +
    geom_vline(xintercept=c(-2,2), color = "grey", size = 1)+ ylim(0,50) +
    geom_text_repel(aes(label = ifelse(abs(log2FoldChange)>1.5 & padj<0.05,geneid,"") ),max.overlaps = Inf, size = 2, segment.size = 0.3 ) +
    ggtitle(sprintf("%s %s vs %s",model,baseline,notbaseline)) + theme_classic()
  
  if (!file.exists(sprintf("%s/Volcano",odir))) {
    dir.create(sprintf("%s/Volcano",odir))
  }
  pdf(sprintf("%s/Volcano/DESeq.pdf", odir))
  print(g)
  graphics.off()
}

## Save RDS object of things
#FIXME: save option
if (TRUE) {
  write.table(cbind(data.frame(contig = rownames(counts_norm)), counts_norm), sprintf("%s/counts_norm.tab", odir),
              sep = "\t", quote = F, row.names = F, col.names = T)
  write.table(outres, sprintf("%s/deseq_res.txt", odir), sep="\t",
              quote=F, row.names = F, col.names = T)
}



