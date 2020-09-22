library(Seurat)
library(ggpubr)
library(ggplot2)

nCoV.integrated <- readRDS(file = "./output/nCoV_scale.rds")
DefaultAssay(nCoV.integrated) <- "integrated"
Idents(nCoV.integrated)

#color
mylevels <- c("macrophages","mDC","B cells", "T cells" , "Doublets" ,"NK cells" ,"Ciliated" ,"neutrophils" ,
              "pDC", "Secretory","Mast cells","plasma cells")
mycolors = c("#9ADFBF","#FFD699","#F7BAA8","#BFBFBF","#F42E3C","#F49FAC","#934CCA","#5E75BA","#90EE90","#3DCFFF","#E2DFA8","#BF8282")
names(mycolors) <- mylevels

nCoV.integrated <- FindClusters(object = nCoV.integrated, resolution = 0.4) 
nCoV.integrated <- FindClusters(object = nCoV.integrated, resolution = 0.8) 
nCoV.integrated <- FindClusters(object = nCoV.integrated, resolution = 1.6) 
nCoV.integrated <- FindClusters(object = nCoV.integrated, resolution = 2) 

p1 <- DimPlot(object = nCoV.integrated, reduction = 'tsne',label = TRUE,group.by = "label", cols = mycolors) + labs(title = "annotation") 
p2 <- DimPlot(object = nCoV.integrated, reduction = 'tsne',label = TRUE,group.by = "integrated_snn_res.0.4")  + labs(title = "resolution=0.4") 
p3 <- DimPlot(object = nCoV.integrated, reduction = 'tsne',label = TRUE,group.by = "integrated_snn_res.0.8")  + labs(title = "resolution=0.8") 
p4 <- DimPlot(object = nCoV.integrated, reduction = 'tsne',label = TRUE,group.by = "integrated_snn_res.1.2")  + labs(title = "resolution=1.2") 
p5 <- DimPlot(object = nCoV.integrated, reduction = 'tsne',label = TRUE,group.by = "integrated_snn_res.1.6")  + labs(title = "resolution=1.6") 
p6 <- DimPlot(object = nCoV.integrated, reduction = 'tsne',label = TRUE,group.by = "integrated_snn_res.2")  + labs(title = "resolution=2") 

p_merge <- ggarrange(p1,p2,p3,p4,p5,p6,ncol = 3,nrow=2,align = "hv")
dpi=500
png(file="resolution/merge.png", width = dpi*25, height = dpi*15, units = "px",res = dpi)
print(p_merge)
dev.off()

################################################################################

#markers
markers <- c("FABP4","CD52","MARCO",   ##"FABP4+ Macro"
               "S100A8","FCN1","S100A9",  ##"FCN1+ Macro"
               "SPP1","LGMN","CCL18",   ##"LGMN+ Macro"
               "CCR7","IL7R",  ##"CCR7+ T"
               "GZMA","GZMB",  ##"CD8+ T"
               "CTLA4","IL2RA",  ##"mixed T"
               "MKI67","TYMS",  ##"Proliferating T"
               'TPPP3','KRT18',  ##"Epi"
               'IL3RA',"LILRA4",  ##"pDC"
               'CD79A','MS4A1',  ##"B"
               'TPSB2',  ##"Mast"
               'IGHG4', ##"Plasma"
               'CD3D',  ##"T"
               'KLRF1',   ##"NK"
               "FCGR3B",  ##"Neutrophils"
               'CD1C','CLEC9A',  ##"mDC"
               'CD14','CD68')  ##"Macro"

#plot
plot_marker <- function(seu.obj,markers,res){
  
  dpi=300
  
  pdf(file=paste("resolution/res",res,"/marker/marker_heatmap.pdf",sep = ""), width = 16, height = 20)
  pp = DotPlot(seu.obj, features = rev(markers),cols = c('white','#F8766D'),dot.scale =5) + RotatedAxis()
  pp = pp + theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12)) + labs(x='',y='') +
    guides(color = guide_colorbar(title = 'Scale expression'),size = guide_legend(title = 'Percent expressed')) +
    theme(axis.line = element_line(size = 0.6))
  print(pp)
  dev.off()
  
  for(marker in markers){
    png(file=paste("resolution/res",res,"/marker/violin_",marker,".png",sep=''), width = dpi*8, height = dpi*3, units = "px",res = dpi)
    print(VlnPlot(object = seu.obj, features = marker,pt.size = 0)+ theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)))
    dev.off()
  }
  
  for(marker in markers){
    png(file=paste("resolution/res",res,"/marker/tsne_",marker,".png",sep=''), width = dpi*6, height = dpi*4, units = "px",res = dpi)
    print(FeaturePlot(object = seu.obj, reduction = "tsne", features = marker))
    dev.off()
  }
  
}

#color
mylevels <- c("Macrophages","mDC","B cells", "T cells" , "Doublets1" ,"NK cells" ,"Ciliated cells" ,
              "Neutrophils","pDC", "Secretory cells","Mast cells","Plasma cells",
              "FABP4+ Macro","FCN1+ Macro","LGMN+ Macro",
              "CCR7+ T","CD8+ T","mixed T","Proliferating T",
              "Doublets2")
mycolors = c("#9ADFBF","#FFD699","#F7BAA8","#BFBFBF","#F42E3C","#F49FAC","#934CCA",
             "#5E75BA","#90EE90","#3DCFFF","#E2DFA8","#BF8282",
             "#B3EE3A","#20B2AA","#CDC673",
             "#8C4545","#E0AEAE","#BB7C7C","#E0D1D1",
             "#C9903B")
names(mycolors) <- mylevels


################################################################################

##resolution=0.4

Idents(nCoV.integrated) <- nCoV.integrated$integrated_snn_res.0.4
DefaultAssay(nCoV.integrated) <- "RNA"

plot_marker(nCoV.integrated,markers,"04")

##anno
cell.ident <- nCoV.integrated$integrated_snn_res.0.4
cell.ident <- as.vector(cell.ident)
cell.ident[cell.ident %in% c(3,10,12)] <- "T cells"
cell.ident[cell.ident %in% c(1,4,5,6)] <- "FABP4+ Macro"
cell.ident[cell.ident %in% c(0,8,15)] <- "FCN1+ Macro"
cell.ident[cell.ident %in% c(2)] <- "LGMN+ Macro"
cell.ident[cell.ident %in% c(14)] <- "Macrophages"
cell.ident[cell.ident %in% c(17)] <- "Doublets1"
cell.ident[cell.ident %in% c(11)] <- "Doublets2"
cell.ident[cell.ident %in% c(9)] <- "Secretory cells"
cell.ident[cell.ident %in% c(7)] <- "Ciliated cells"
cell.ident[cell.ident %in% c(13)] <- "Plasma cells"
cell.ident[cell.ident %in% c(18,19)] <- "pDC"
cell.ident[cell.ident %in% c(16)] <- "Neutrophils"
names(cell.ident) <- names(nCoV.integrated$integrated_snn_res.0.4)
table(cell.ident)
nCoV.integrated$label04 <- cell.ident


p7 <- DimPlot(object = nCoV.integrated, reduction = 'tsne',label = TRUE) + 
  labs(title = "resolution=0.4")  +
  theme(plot.title = element_text(hjust = 0.5))
p8 <- DimPlot(object = nCoV.integrated, reduction = 'tsne',label = TRUE,group.by = "label04",cols = mycolors,repel = T,label.size = 3) + 
  labs(title = "resolution=0.4")  +
  theme(plot.title = element_text(hjust = 0.5))
p78 <- ggarrange(p7,p8,ncol = 2,align = "hv")

dpi=500
png(file="resolution/res04/res04.png", width = dpi*15, height = dpi*6, units = "px",res = dpi)
print(p78)
dev.off()

barcodetype <- data.frame(Barcode=colnames(nCoV.integrated)[nCoV.integrated$disease == "Y"],
                          Cluster=nCoV.integrated$label04[nCoV.integrated$disease == "Y"],
                          stringsAsFactors = F)
write.table(barcodetype,"resolution/res04/res04_barcodetype.txt",quote = F,sep = "\t",row.names = F)

################################################################################

##resolution=0.8

Idents(nCoV.integrated) <- nCoV.integrated$integrated_snn_res.0.8
DefaultAssay(nCoV.integrated) <- "RNA"

plot_marker(nCoV.integrated,markers,"08")

##anno
cell.ident <- nCoV.integrated$integrated_snn_res.0.8
cell.ident <- as.vector(cell.ident)
cell.ident[cell.ident %in% c(7,12,18)] <- "T cells"
cell.ident[cell.ident %in% c(9)] <- "CD8+ T"
cell.ident[cell.ident %in% c(1,3,4,8)] <- "FABP4+ Macro"
cell.ident[cell.ident %in% c(0,2,10,20,21,22)] <- "FCN1+ Macro"
cell.ident[cell.ident %in% c(5,6)] <- "LGMN+ Macro"
cell.ident[cell.ident %in% c(16)] <- "Doublets2"
cell.ident[cell.ident %in% c(17)] <- "Doublets1"
cell.ident[cell.ident %in% c(11)] <- "Secretory cells"
cell.ident[cell.ident %in% c(13,15)] <- "Ciliated cells"
cell.ident[cell.ident %in% c(19,24)] <- "Plasma cells"
cell.ident[cell.ident %in% c(23)] <- "pDC"
cell.ident[cell.ident %in% c(14)] <- "Neutrophils"
names(cell.ident) <- names(nCoV.integrated$integrated_snn_res.0.8)
table(cell.ident)
nCoV.integrated$label08 <- cell.ident


p9 <- DimPlot(object = nCoV.integrated, reduction = 'tsne',label = TRUE) + 
  labs(title = "resolution=0.8")  +
  theme(plot.title = element_text(hjust = 0.5))
p10 <- DimPlot(object = nCoV.integrated, reduction = 'tsne',label = TRUE, group.by = "label08",cols = mycolors,repel = T,label.size = 3) + 
  labs(title = "resolution=0.8")  +
  theme(plot.title = element_text(hjust = 0.5))
p910 <- ggarrange(p9,p10,ncol = 2,align = "hv")

dpi=500
png(file="resolution/res08/res08.png", width = dpi*15, height = dpi*6, units = "px",res = dpi)
print(p910)
dev.off()

barcodetype <- data.frame(Barcode=colnames(nCoV.integrated)[nCoV.integrated$disease == "Y"],
                          Cluster=nCoV.integrated$label08[nCoV.integrated$disease == "Y"],
                          stringsAsFactors = F)
write.table(barcodetype,"resolution/res08/res08_barcodetype.txt",quote = F,sep = "\t",row.names = F)

################################################################################

##resolution=1.2

Idents(nCoV.integrated) <- nCoV.integrated$integrated_snn_res.1.2
DefaultAssay(nCoV.integrated) <- "RNA"

plot_marker(nCoV.integrated,markers,"12")

##anno
cell.ident <- nCoV.integrated$integrated_snn_res.1.2
cell.ident <- as.vector(cell.ident)
cell.ident[cell.ident %in% c(8,15)] <- "T cells"
cell.ident[cell.ident %in% c(10)] <- "CD8+ T"
cell.ident[cell.ident %in% c(3,4,9,11)] <- "FABP4+ Macro"
cell.ident[cell.ident %in% c(0,1,2,12,16,22,23,6)] <- "FCN1+ Macro"
cell.ident[cell.ident %in% c(5,7,19)] <- "LGMN+ Macro"
cell.ident[cell.ident %in% c(21)] <- "Doublets1"
cell.ident[cell.ident %in% c(14)] <- "Secretory cells"
cell.ident[cell.ident %in% c(18,26,29)] <- "Ciliated cells"
cell.ident[cell.ident %in% c(24,25,31)] <- "Plasma cells"
cell.ident[cell.ident %in% c(28)] <- "pDC"
cell.ident[cell.ident %in% c(13)] <- "Neutrophils"
cell.ident[cell.ident %in% c(17)] <- "NK cells"
cell.ident[cell.ident %in% c(20)] <- "mDC"
cell.ident[cell.ident %in% c(27)] <- "B cells"
cell.ident[cell.ident %in% c(30)] <- "Mast cells"
names(cell.ident) <- names(nCoV.integrated$integrated_snn_res.1.2)
table(cell.ident)
nCoV.integrated$label12 <- cell.ident


p11 <- DimPlot(object = nCoV.integrated, reduction = 'tsne',label = TRUE) + 
  labs(title = "resolution=1.2")  +
  theme(plot.title = element_text(hjust = 0.5))
p12 <- DimPlot(object = nCoV.integrated, reduction = 'tsne',label = TRUE, group.by = "label12",cols = mycolors,repel = T,label.size = 3) + 
  labs(title = "resolution=1.2")  +
  theme(plot.title = element_text(hjust = 0.5))
p1112 <- ggarrange(p10,p11,ncol = 2,align = "hv")

dpi=500
png(file="resolution/res12/res12.png", width = dpi*15, height = dpi*6, units = "px",res = dpi)
print(p1112)
dev.off()

barcodetype <- data.frame(Barcode=colnames(nCoV.integrated)[nCoV.integrated$disease == "Y"],
                          Cluster=nCoV.integrated$label12[nCoV.integrated$disease == "Y"],
                          stringsAsFactors = F)
write.table(barcodetype,"resolution/res12/res12_barcodetype.txt",quote = F,sep = "\t",row.names = F)

################################################################################

##resolution=1.6

Idents(nCoV.integrated) <- nCoV.integrated$integrated_snn_res.1.6
DefaultAssay(nCoV.integrated) <- "RNA"

plot_marker(nCoV.integrated,markers,"16")

##anno
cell.ident <- nCoV.integrated$integrated_snn_res.1.6
cell.ident <- as.vector(cell.ident)
cell.ident[cell.ident %in% c(10)] <- "CCR7+ T"
cell.ident[cell.ident %in% c(11)] <- "CD8+ T"
cell.ident[cell.ident %in% c(16)] <- "Proliferating T"
cell.ident[cell.ident %in% c(1,13,17,19,27,3,4,9)] <- "FABP4+ Macro"
cell.ident[cell.ident %in% c(0,12,14,23,24,5,7)] <- "FCN1+ Macro"
cell.ident[cell.ident %in% c(2,21,6,8)] <- "LGMN+ Macro"
cell.ident[cell.ident %in% c(25)] <- "Doublets1"
cell.ident[cell.ident %in% c(15)] <- "Secretory cells"
cell.ident[cell.ident %in% c(18,29,32)] <- "Ciliated cells"
cell.ident[cell.ident %in% c(26,28,34)] <- "Plasma cells"
cell.ident[cell.ident %in% c(31)] <- "pDC"
cell.ident[cell.ident %in% c(14)] <- "Neutrophils"
cell.ident[cell.ident %in% c(20)] <- "NK cells"
cell.ident[cell.ident %in% c(22)] <- "mDC"
cell.ident[cell.ident %in% c(30)] <- "B cells"
cell.ident[cell.ident %in% c(33)] <- "Mast cells"
names(cell.ident) <- names(nCoV.integrated$integrated_snn_res.1.6)
table(cell.ident)
nCoV.integrated$label16 <- cell.ident


p13 <- DimPlot(object = nCoV.integrated, reduction = 'tsne',label = TRUE)+ 
  labs(title = "resolution=1.6")  +
  theme(plot.title = element_text(hjust = 0.5)) 
p14 <- DimPlot(object = nCoV.integrated, reduction = 'tsne',label = TRUE, group.by = "label16",cols = mycolors,repel = T,label.size = 3) + 
  labs(title = "resolution=1.6")  +
  theme(plot.title = element_text(hjust = 0.5))
p1314 <- ggarrange(p12,p13,ncol = 2,align = "hv")

dpi=500
png(file="resolution/res16/res16.png", width = dpi*15, height = dpi*6, units = "px",res = dpi)
print(p1314)
dev.off()

barcodetype <- data.frame(Barcode=colnames(nCoV.integrated)[nCoV.integrated$disease == "Y"],
                          Cluster=nCoV.integrated$label16[nCoV.integrated$disease == "Y"],
                          stringsAsFactors = F)
write.table(barcodetype,"resolution/res16/res16_barcodetype.txt",quote = F,sep = "\t",row.names = F)


################################################################################

##resolution=2

Idents(nCoV.integrated) <- nCoV.integrated$integrated_snn_res.2
DefaultAssay(nCoV.integrated) <- "RNA"

plot_marker(nCoV.integrated,markers,"2")

##anno
cell.ident <- nCoV.integrated$integrated_snn_res.2
cell.ident <- as.vector(cell.ident)
cell.ident[cell.ident %in% c(31)] <- "mixed T"
cell.ident[cell.ident %in% c(7)] <- "CCR7+ T"
cell.ident[cell.ident %in% c(8)] <- "CD8+ T"
cell.ident[cell.ident %in% c(17)] <- "Proliferating T"
cell.ident[cell.ident %in% c(0,2,3,9,15,20,21,22,26,28,32,37)] <- "FABP4+ Macro"
cell.ident[cell.ident %in% c(4,5,10,11,12,14,29)] <- "FCN1+ Macro"
cell.ident[cell.ident %in% c(1,6,13,18,19)] <- "LGMN+ Macro"
cell.ident[cell.ident %in% c(25)] <- "Doublets1"
cell.ident[cell.ident %in% c(16)] <- "Secretory cells"
cell.ident[cell.ident %in% c(24,35,39)] <- "Ciliated cells"
cell.ident[cell.ident %in% c(33,34,41)] <- "Plasma cells"
cell.ident[cell.ident %in% c(38)] <- "pDC"
cell.ident[cell.ident %in% c(30)] <- "Neutrophils"
cell.ident[cell.ident %in% c(23)] <- "NK cells"
cell.ident[cell.ident %in% c(27)] <- "mDC"
cell.ident[cell.ident %in% c(36)] <- "B cells"
cell.ident[cell.ident %in% c(40)] <- "Mast cells"
names(cell.ident) <- names(nCoV.integrated$integrated_snn_res.2)
table(cell.ident)
nCoV.integrated$label2 <- cell.ident


p15 <- DimPlot(object = nCoV.integrated, reduction = 'tsne',label = TRUE) + 
  labs(title = "resolution=2")  +
  theme(plot.title = element_text(hjust = 0.5))
p16 <- DimPlot(object = nCoV.integrated, reduction = 'tsne',label = TRUE, group.by = "label2",cols = mycolors,repel = T,label.size = 3) + 
  labs(title = "resolution=2")  +
  theme(plot.title = element_text(hjust = 0.5))
p1516 <- ggarrange(p14,p15,ncol = 2,align = "hv")

dpi=500
png(file="resolution/res2/res2.png", width = dpi*15, height = dpi*6, units = "px",res = dpi)
print(p1516)
dev.off()

barcodetype <- data.frame(Barcode=colnames(nCoV.integrated)[nCoV.integrated$disease == "Y"],
                          Cluster=nCoV.integrated$label2[nCoV.integrated$disease == "Y"],
                          stringsAsFactors = F)
write.table(barcodetype,"resolution/res2/res2_barcodetype.txt",quote = F,sep = "\t",row.names = F)

#################################################################################################

p_merge <- ggarrange(p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,
                     ncol=2,nrow = 5,align = "hv",font.label = list(size=24),
                     labels = c("A","","B","","C","","D","","E",""))

dpi=500
tiff(file="resolution/merge_tsne.tiff", width = dpi*12, height = dpi*21.5, units = "px",res = dpi)
print(p_merge)
dev.off()

saveRDS(nCoV.integrated,"./resolution/nCoV_res.rds")
