library(tidyverse)
library(ggrepel)

#### Load in data
RNA<-read_tsv(snakemake@input[["tsv_data"]])
RNA<-RNA%>%filter(baseMean>=snakemake@input[["baseline"]])%>%na.omit()

#### Metasheet
metasheet<-function(up,down){
  no_change<-nrow(RNA[between(log2FoldChange,down,up)])
}


no_change<-RNA%>%filter(between(log2FoldChange,-1,1))%>%nrow()
increase<-RNA%>%filter(log2FoldChange>=1)%>%filter(pvalue<=0.05)%>%nrow()
decrease<-RNA%>%filter(log2FoldChange<=-1)%>%filter(pvalue<=0.05)%>%nrow()

stats<-data.frame(
  no_change = no_change, 
  increase = increase,
  decrease = decrease)

write.table(stats,file=paste0(snakemake@params[["odir"]],"stats.txt"),sep='\t')

#### Top genes

top_up<-head(RNA[order(RNA$log2FoldChange,decreasing = T),],
             snakemake@params[["top_number"]])
top_down<-head(RNA[order(RNA$log2FoldChange,decreasing = F),],
               snakemake@params[["top_number"]])
wgt_top_up<-head(RNA[order(RNA$log2FoldChange*RNA$baseMean,decreasing = T),],
                 snakemake@params[["top_number"]])
wgt_top_down<-head(RNA[order(RNA$log2FoldChange*RNA$baseMean,decreasing = F),],
                   snakemake@params[["top_number"]])

top_genes<-rbind(top_ten_up,top_ten_down,wgt_top_ten_up,wgt_top_ten_down)

write.table(stats,file=paste0(snakemake@params[["odir"]],"top_genes.txt"),sep='\t')


#### Volcano Plot
label_pos<-seq(50,25,-1.5) ## sets y coordinate for annotation


plot<-ggplot()+
  geom_point(data=RNA,
             aes(x=log2FoldChange,y=-log10(padj)),alpha=0.6,
             color=case_when(RNA$log2FoldChange< -1 & RNA$padj<0.05~"#B22222",
                             RNA$log2FoldChange>1&RNA$padj<0.05~"#004D40",
                             T~"black"))+
  geom_hline(yintercept = -log10(0.05),color="#4682B4",size=1,linetype=2)+
  xlab(snakemake@input[["contrast"]])+ ylab("-log10 padj")+
  theme_classic()+
  coord_cartesian(xlim=c(-6,6),ylim=c(0,50))+
  theme(axis.text.x = element_text(size =20,color="black"),
        axis.text.y = element_text(size =20,color="black"),
        axis.title.x = element_text(size =18,color="black"),
        axis.title.y = element_text(size =18,color="black"),
        axis.line = element_line(linewidth=1.5),axis.ticks = element_line(linewidth=1.5),
        axis.ticks.length = unit(0.2,"cm"))+
  annotate("text", x=-4,y=label_pos[1:11],label=c("Top Decreased Genes\n",top_down$hgnc_symbol[1:10]),
           hjust = 0)+
  annotate("text", x=2,y=label_pos[1:11],label=c("Top Increased Genes\n",top_up$hgnc_symbol[1:10]),
           hjust = 0)


ggsave(plot, "svg", path=paste0(snakemake@output[["plot"]],"volcano.svg"), dpi="print",width=7,height=5, units = "in")
