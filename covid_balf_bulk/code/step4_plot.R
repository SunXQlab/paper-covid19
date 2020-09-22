options(stringsAsFactors = F)
library(ggplot2)

##TF
TF_res=read.table("output/fisher_TF.txt",sep = "\t")
TF_res$neglog10_padj <- -log10(TF_res$pval)
TF_res$type <- rep("TF",nrow(TF_res))

p1 <- ggplot(data=TF_res, aes(x=neglog10_padj, y=reorder(TF,neglog10_padj),fill=type)) +
  geom_bar(stat="identity", width=0.8) + xlim(c(0,10)) +
  scale_fill_manual(values = "#FD8D62") + theme_test()  +
  labs(title = "The activate Transcription  Factor",x="-log10(pval)",y="") + 
  theme(axis.text=element_text(face = "bold", color="gray50"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label=reorder(TF,neglog10_padj), x=neglog10_padj), hjust = -0.1,vjust=0.5) 

##pathway
pw_res=read.table("output/fisher_pathway.txt",sep = "\t")
pw_res$neglog10_padj <- -log10(pw_res$pval)
pw_res$type <- rep("kegg",nrow(pw_res))
pw_res <- pw_res[order(pw_res$neglog10_padj),]
pw_res$pathway <- factor(pw_res$pathway,levels = pw_res$pathway)

p2 <- ggplot(data=pw_res, aes(x=neglog10_padj, y=pathway,fill=type)) +
  geom_bar(stat="identity", width=0.8) + xlim(c(0,12)) +
  scale_fill_manual(values = "#FD8D62") + theme_test()  +
  labs(title = "The activate KEGG pathway",x="-log10(pval)",y="") + 
  theme(axis.text=element_text(face = "bold", color="gray50"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label=pathway, x=neglog10_padj), hjust = -0.1,vjust=0.5) 

##merge
p3 <- ggarrange(p1,p2, ncol = 2, labels = c("A", "B"),align = "hv",common.legend = F)
ggsave(p3,filename = "fisher_TFPW.pdf",width = 12, height = 8)

########################################################################################

#heatmap
library(pheatmap)

expr_norm <- read.table(file='data/norm.txt',header = T,sep='\t')
group <- factor(c(rep("Normal",3),rep("Patient",4)))
group <- relevel(group, ref = "Normal")

deg=read.table('./output/DE_genes.txt',sep = '\t',header = T)
deg=na.omit(deg)
deg_sig=deg[abs(deg$log2FoldChange) > log2(4) & deg$padj < 0.01,]
deg_sig=deg_sig[order(deg_sig$log2FoldChange),]

deg_up=head(deg_sig$gene,100)
deg_down=tail(deg_sig$gene,100)

hvg=c(deg_up,deg_down)
hvg_scale=t(scale(t(expr_norm[hvg,]))) 
hvg_scale[hvg_scale>2]=2 
hvg_scale[hvg_scale< -2]= -2

group <- factor(c(rep("Normal",3),rep("Patient",4)))
ac=data.frame(g=group)
rownames(ac)=colnames(hvg_scale)

pdf("pheatmap_path_gene.pdf",width = 6,height = 8)
pheatmap(hvg_scale,show_colnames =T,show_rownames = F,
         cluster_cols = F,cluster_rows = T,
         treeheight_row = 0,treeheight_col = 0,
         annotation_col=ac,angle_col = "0",
         border_color = NA,fontsize_row = 8)
dev.off()

