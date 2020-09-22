##integrate data with seurat v3
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)

samples = read.delim2("data/meta.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
nCoV.list = list()
for(sample_s in samples$sample){
  print(sample_s)
  sample_i = samples %>% dplyr::filter(.,sample == sample_s)
  
  if(sample_s == "GSM3660650"){
    datadir = paste0("data/scRNA-seq/",sample_s,"/filtered_feature_bc_matrix")
    nCoV.data.i <- Read10X(data.dir = datadir)
  }else{
    datadir = paste0("data/scRNA-seq/",sample_s,"_filtered_feature_bc_matrix.h5")
    nCoV.data.i <- Read10X_h5(filename = datadir)
  }
  
  sample.tmp.seurat <- CreateSeuratObject(counts = sample.tmp, min.cells = 3, min.features = 200,project = sample_s)
  sample.tmp.seurat[['percent.mito']] <- PercentageFeatureSet(sample.tmp.seurat, pattern = "^MT-")
  sample_i$nFeature_RNA_low = as.numeric(sample_i$nFeature_RNA_low)
  sample_i$nFeature_RNA_high = as.numeric(sample_i$nFeature_RNA_high)
  sample_i$nCount_RNA = as.numeric(sample_i$nCount_RNA)
  sample_i$percent.mito = as.numeric(sample_i$percent.mito)
  sample.tmp.seurat <- subset(x = sample.tmp.seurat, subset = nFeature_RNA > sample_i$nFeature_RNA_low & nFeature_RNA < sample_i$nFeature_RNA_high 
                              & nCount_RNA > sample_i$nCount_RNA & percent.mito < sample_i$percent.mito)
  sample.tmp.seurat <- NormalizeData(sample.tmp.seurat, verbose = FALSE)
  sample.tmp.seurat <- FindVariableFeatures(sample.tmp.seurat, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
  nCoV.list[sample_s] = sample.tmp.seurat
}
nCoV <- FindIntegrationAnchors(object.list = nCoV.list, dims = 1:50)
nCoV.integrated <- IntegrateData(anchorset = nCoV, dims = 1:50,features.to.integrate = rownames(nCoV))

####add  sample info
sample_info = as.data.frame(colnames(nCoV.integrated))
colnames(sample_info) = c('ID')
rownames(sample_info) = sample_info$ID
sample_info$sample = nCoV.integrated@meta.data$orig.ident
sample_info = dplyr::left_join(sample_info,samples)
rownames(sample_info) = sample_info$ID
nCoV.integrated = AddMetaData(object = nCoV.integrated, metadata = sample_info)

###first generate data and scale data in RNA assay
DefaultAssay(nCoV.integrated) <- "RNA"
nCoV.integrated[['percent.mito']] <- PercentageFeatureSet(nCoV.integrated, pattern = "^MT-")
nCoV.integrated <- NormalizeData(object = nCoV.integrated, normalization.method = "LogNormalize", scale.factor = 1e4)
nCoV.integrated <- FindVariableFeatures(object = nCoV.integrated, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
nCoV.integrated <- ScaleData(nCoV.integrated, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))

##change to integrated assay
DefaultAssay(nCoV.integrated) <- "integrated"
dpi = 300
png(file="output/figure/qc.png", width = dpi*16, height = dpi*8, units = "px",res = dpi)
VlnPlot(object = nCoV.integrated, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
dev.off()

png(file="output/figure/umi-gene.png", width = dpi*6, height = dpi*5, units = "px",res = dpi)
FeatureScatter(object = nCoV.integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

# Run the standard workflow for visualization and clustering
nCoV.integrated <- ScaleData(nCoV.integrated, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))
nCoV.integrated <- RunPCA(nCoV.integrated, verbose = FALSE,npcs = 100)
nCoV.integrated <- ProjectDim(object = nCoV.integrated)
png(file="pca.png", width = dpi*10, height = dpi*6, units = "px",res = dpi)
ElbowPlot(object = nCoV.integrated,ndims = 100)
dev.off()

###cluster
nCoV.integrated <- FindNeighbors(object = nCoV.integrated, dims = 1:50)
nCoV.integrated <- FindClusters(object = nCoV.integrated, resolution = 1.2) 

###tsne and umap
nCoV.integrated <- RunTSNE(object = nCoV.integrated, dims = 1:50)
nCoV.integrated <- RunUMAP(nCoV.integrated, reduction = "pca", dims = 1:50)
png(file="output/figure/tsne.png", width = dpi*8, height = dpi*6, units = "px",res = dpi)
DimPlot(object = nCoV.integrated, reduction = 'tsne',label = TRUE)
dev.off()
png(file="output/figure/umap.png", width = dpi*8, height = dpi*6, units = "px",res = dpi)
DimPlot(object = nCoV.integrated, reduction = 'umap',label = TRUE)
dev.off()

dpi = 300
png(file="output/figure/feature.png", width = dpi*24, height = dpi*5, units = "px",res = dpi)
VlnPlot(object = nCoV.integrated, features = c("nFeature_RNA", "nCount_RNA"))
dev.off()

#save seurat object
saveRDS(nCoV.integrated, file = "output/nCoV.rds")

########################################################################################################

library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(reshape2)

##Markers
DefaultAssay(nCoV.integrated) <- "RNA"
# find markers for every cluster compared to all remaining cells, report only the positive ones
nCoV.integrated@misc$markers <- FindAllMarkers(object = nCoV.integrated, assay = 'RNA',only.pos = TRUE, test.use = 'MAST')
write.table(nCoV.integrated@misc$markers,file='output/marker/marker_MAST.txt',row.names = FALSE,quote = FALSE,sep = '\t')

##heatmap
hc.markers = read.delim2("output/marker/marker_MAST.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
hc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10
tt1 = DoHeatmap(object = subset(nCoV.integrated, downsample = 500), features = top10$gene) + NoLegend()
ggplot2::ggsave(file="output/marker/marker_heatmap_MAST.pdf",plot = tt1,device = 'pdf',width = 20, height = 16, units = "in",dpi = dpi,limitsize = FALSE)

##re-scale data in RNA assay and save seurat object
hc.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC) -> top30
var.genes = c(nCoV.integrated@assays$RNA@var.features,top30$gene,markers)

nCoV.integrated <- ScaleData(nCoV.integrated, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"),features = var.genes)
saveRDS(nCoV.integrated, file = "output/nCoV_scale.rds")

####marker expression
dpi = 300
markers = c('AGER','SFTPC','SCGB3A2','TPPP3','KRT5',
            'CD68','FCN1','CD1C','TPSB2','CD14','MARCO','CXCR2',
            'CLEC9A','IL3RA',
            'CD3D','CD8A','KLRF1',
            'CD79A','IGHG4','MS4A1',
            'VWF','DCN',
            'FCGR3A','TREM2','KRT18','HBB')

for(marker in markers){
  png(file=paste("output/marker/violin_",marker,".png",sep=''), width = dpi*8, height = dpi*3, units = "px",res = dpi)
  print(VlnPlot(object = nCoV.integrated, features = marker,pt.size = 0)+ theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)))
  dev.off()

  png(file=paste("output/marker/umap_",marker,".png",sep=''), width = dpi*6, height = dpi*4, units = "px",res = dpi)
  print(FeaturePlot(object = nCoV.integrated, features = marker,cols = c("lightgrey","#ff0000")))
  dev.off()
}

library(ggplot2)
pdf(file="output/marker_heatmap.pdf", width = 10, height = 8)
pp = DotPlot(nCoV.integrated, features = rev(markers),cols = c('white','#F8766D'),dot.scale =5) + RotatedAxis()
pp = pp + theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12)) + labs(x='',y='') +
  guides(color = guide_colorbar(title = 'Scale expression'),size = guide_legend(title = 'Percent expressed')) +
  theme(axis.line = element_line(size = 0.6))
print(pp)
dev.off()