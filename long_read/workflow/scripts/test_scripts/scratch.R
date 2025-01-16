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
