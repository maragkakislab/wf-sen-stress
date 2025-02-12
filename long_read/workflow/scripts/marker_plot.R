library(tidyverse)
library(ggrepel)

#### Load in data
RNA<-read_tsv(snakemake@input[["tsv_data"]])
markers<-snakemake@params[["marker_genes"]]

print(markers)

df<-RNA%>%filter(baseMean>=snakemake@params[["baseline"]])%>%na.omit()%>%
filter(hgnc_symbol %in% markers)


####
write.table(markers,file=paste0(snakemake@params[["odir"]],"markers.txt"),sep='\t', row.names=FALSE)

### plot function
sen_markers<-function(df){
  df_sen<-df%>%filter(hgnc_symbol %in% smark)
  plot<-df_sen%>%ggplot()+
    geom_point(aes(x=hgnc_symbol, y = log2FoldChange,
                   color= padj),size=4)+
    scale_color_gradient(low="red",high="blue",limits=c(0,0.05))+
    geom_hline(yintercept = 0)+
    theme_classic()+theme(axis.text.x = element_text(angle=45, vjust = 0.5))+
    theme(axis.text = element_text(size=14,color="black"))+
    coord_cartesian(ylim=c(-10,10))
    
  return(plot)
}

####

ggsave(filename= snakemake@output[["plot"]], plot = sen_markers(df), device = "pdf", width=7,height=5, units = "in")
