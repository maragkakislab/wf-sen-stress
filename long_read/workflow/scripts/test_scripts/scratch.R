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


write.table(comparison,file=paste0(odir,"contrasts.txt"), sep = "\t")

odir = "/data/Maragkakislab/payeamj/sen_stress/"
#####

for (i in resultsNames(dds)) {
  assign(i, as.data.frame(results(dds, name = i)))
}

####


as.data.frame(results(ddsMat, contrast=c("Condition","24H","0H")))

sample<- c("A","B","C","D")
value<-c(1,9,0,1)

df1<-data.frame(sample = sample, value =  value)
df1

df2<-data.frame(t(value))
colnames(df2)<-sample

df2

if (colnames(df2 != data_file$sample)) {
  stop("Count data library names are not in metadata!")
}


ifelse(colnames(df2) == df1$sample,"True","False")

identical(colnames(df2),df1$sample)



z<-as.data.frame((combn(sample,2)))

colnames(z)<-rep("sample", dim(z)[2])

z

for (i in colnames(z)){
  print(paste("sample", z$[[i]][1], z$[[i]][2]))}
for (i in colnames(z)) {
  print(paste("sample", z[[i]][1], z[[i]][2]))
}
