#### pratice deseq2 dataset and metadat from bioconductor page
library("DESeq2")
library("pasilla")
pasCts <- system.file("extdata",
                      "pasilla_gene_counts.tsv",
                      package="pasilla", mustWork=TRUE)
pasAnno <- system.file("extdata",
                       "pasilla_sample_annotation.csv",
                       package="pasilla", mustWork=TRUE)
cts <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))
coldata <- read.csv(pasAnno, row.names=1)
coldata <- coldata[,c("condition","type")]
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)
rownames(coldata) <- sub("fb", "", rownames(coldata))
cts <- cts[, rownames(coldata)]


######

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds <- DESeq(dds)
comparison<-resultsNames(dds)

#######
library(gtools)
metafile <- read.delim("/data/Maragkakislab/payeamj/sen_stress/long_read/analysis/ars_viability/metafile.txt", 
                       header = TRUE, sep = "\t", 
                       stringsAsFactors = F)

conditions<-unique(metafile$condition)
permutations <- permutations(n = length(conditions), 
                             r = 2, 
                             v = conditions)
permutations_list <- split(permutations, row(permutations))

