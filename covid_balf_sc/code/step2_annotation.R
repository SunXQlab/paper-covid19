library(Seurat)
library(ggplot2)
library(ggpubr)
library(reshape2)

##Annotation
#read seurat object
nCoV.integrated <- readRDS(file = "./output/nCoV_scale.rds")

#anno
cell.ident <- nCoV.integrated@active.ident
label <- c(rep("macrophages",8),"T cells","macrophages","T cells","macrophages","macrophages",
           "neutrophils","Secretory","T cells","macrophages","NK cells","Ciliated","macrophages",
           "mDC","Doublets",rep("macrophages",2),rep("plasma cells",2),"Ciliated","B cells",
           "pDC","Ciliated","Mast cells","macrophages")
levels(cell.ident) <- label
nCoV.integrated$label <- cell.ident

mylevels <- c("macrophages","mDC","B cells", "T cells" , "Doublets" ,"NK cells" ,"Ciliated" ,"neutrophils" ,
              "pDC", "Secretory","Mast cells","plasma cells")
mycolors = c("#9ADFBF","#FFD699","#F7BAA8","#BFBFBF","#F42E3C","#F49FAC","#934CCA","#5E75BA","#90EE90","#3DCFFF","#E2DFA8","#BF8282")
names(mycolors) <- mylevels

saveRDS(nCoV.integrated,file = "./output/nCoV_scale.rds")

###############################################################################################

#check anno
orig.anno <- read.table("data/all.cell.annotation.meta.txt",header = T,sep = "\t")

new.anno <- data.frame(ID=names(cell.ident),
                       anno=as.character(cell.ident))
new.anno$ID <- gsub(pattern = "\\W\\d+",replacement = "",x = new.anno$ID)

df_merge <- merge(orig.anno,new.anno,by.x = "ID",all.x = T)
com <- table(df_merge$anno,df_merge$celltype)

df_com <- as.data.frame(com)
df_com <- acast(df_com,Var1~Var2)
write.table(df_com,"./output/compare.txt",quote = F,sep = "\t",row.names = T,col.names = T)

###############################################################################################

#select patient
nCoV.pat <- subset(nCoV.integrated,disease == "Y")

#marker
markers = c("CD68", #macro
            "TPSB2", #Mast cells
            "FCGR3B","FCGR3A", #neutrophils
            "CLEC9A","CD1C", #mDCs
            "LILRA4","IL3RA", #pDCs
            "KLRD1", #NK
            "CD3D", #T cells
            "MS4A1", #B cells
            "IGHG4", # plasma cells
            "TPPP3", "KRT18") # epithelial cells

#pdf(file="./figure/marker_heatmap.pdf", width = 10, height = 4)
pp = DotPlot(nCoV.pat, features = rev(markers),cols = c('white','#F8766D'),group.by = "label",dot.scale =5) + RotatedAxis()
pp = pp + theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12)) + labs(x='',y='') +
  guides(color = guide_colorbar(title = 'Scale expression'),size = guide_legend(title = 'Percent expressed')) +
  theme(axis.line = element_line(size = 0.6))
print(pp)
#dev.off()

dpi=300
for(marker in markers){
  png(file=paste("figure/marker/violin_",marker,".png",sep=''), width = dpi*8, height = dpi*3, units = "px",res = dpi)
  print(VlnPlot(object = nCoV.pat, features = marker,pt.size = 0, group.by = "label")+ theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)))
  dev.off()
}

plot_ls <- list()
for(i in 1:length(markers)){
  
  marker <- markers[i]
  plot_ls[[marker]] <- FeaturePlot(object = nCoV.pat, reduction = "tsne", features = marker)
  
}

p_merge<-ggarrange(plotlist = plot_ls,
          nrow = 4,ncol = 4,align = "hv",
          font.label = list(size=16),common.legend = F)

png(file=paste("figure/tsne_marker.png",sep=''), width = dpi*24, height = dpi*24, units = "px",res = dpi)
print(p_merge)
dev.off()

tiff(file=paste("figure/tsne_marker.tiff",sep=''),width = 2400, height = 2400)
print(p_merge)
dev.off()

################################################################################
#SARA-CoV2 receptor

#run ALRA
nCoV.alra <- RunALRA(nCoV.integrated)
#save seurat object
saveRDS(nCoV.alra,file = "./output/nCoV_alra.rds")

nCoV.pat <- subset(nCoV.alra,disease == "Y")

p1 <- DimPlot(object = nCoV.pat, reduction = 'tsne',label = FALSE, group.by = "label",cols = mycolors)+
  theme(legend.position = "bottom") 

p2 <- FeaturePlot(object = nCoV.alra, reduction = 'tsne', features = "ACE2",cols = c("lightgrey","#ff0000")) +
  labs(colour = "ACE2 expression") +
  theme(plot.title = element_blank(),
        legend.text.align = 0,
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.justification = c(0.5,0.5)) 

p3 <-VlnPlot(object = nCoV.alra,features = "ACE2",pt.size = 0, group.by = "label",cols = mycolors) +
  xlab("Cell Type") + ylab("ACE2 Expression") +
  theme(axis.text.x = element_text(size = 12,angle = 0,hjust = 0.5),
        legend.position = "none",plot.title = element_blank())

p4 <- ggarrange(p1,p2,ncol = 2,align = "hv") 
p5 <- ggarrange(p4,p3,nrow = 2,heights = c(2,1))

dpi=500
png(file="figure/merge.png", width = dpi*12, height = dpi*10.5, units = "px",res = dpi)
print(p5)
dev.off()

pdf(file="figure/merge.pdf", width = 12, height = 10.5)
print(p5)
dev.off()

################################################################################

p6 <- DimPlot(object = nCoV.pat, reduction = 'tsne',label = FALSE, group.by = "group")

p7 <- FeaturePlot(object = nCoV.alra, reduction = 'tsne', features = "ACE2",split.by = "group") +
  theme(legend.position = "none") 

p8 <- ggarrange(ggarrange(pp,p7,nrow = 2,align = "hv",common.legend = F,labels = c("A","B"))
                ,p6,
                nrow = 2,heights = c(2,1),align = "hv",common.legend = F,labels = c("","C"))

pdf(file="figure/scrnaseq_detail.pdf", width = 12, height = 10.5)
print(p8)
dev.off()

################################################################################

#prepare the input for scMLnet

#patient data
pat.barcodetype <- barcodetype[barcodetype$Barcode %in% colnames(nCoV.pat),]
write.table(pat.barcodetype,"./data/pat.barcodetype.txt",quote = F,sep = "\t",row.names = F)

pat.rawcount <- nCoV.pat@assays$RNA@counts
saveRDS(pat.rawcount,file = "./data/pat_data.Rdata")
