#function
getCorR <- function(GCMat,keys,cells)
{
  
  GCMat <- GCMat[,cells]
  GCMat <- as.matrix(GCMat)
  
  #check key
  keys <- unique(keys)
  key_df <- data.frame(key=keys)
  key_df$gene.1 <- unlist(lapply(keys,function(key){strsplit(key,"_")[[1]][1]}))
  key_df$gene.2 <- unlist(lapply(keys,function(key){strsplit(key,"_")[[1]][2]}))
  
  key_df <- key_df[ key_df$gene.1 %in% rownames(GCMat) & key_df$gene.2 %in% rownames(GCMat),]
  key_df <- key_df[rowMeans(GCMat[key_df$gene.1,])>0 & rowMeans(GCMat[key_df$gene.2,])>0 ,]
  
  #cal cor
  cor_df <- matrix(nrow = nrow(key_df), ncol = 2)
  for(i in 1:nrow(key_df)){
    
    key_info <- key_df[i,]
    
    cor <- cor.test(GCMat[key_info$gene.1,],GCMat[key_info$gene.2,],method = c("kendall"))
    pval <- signif(cor$p.value,digits = 3)
    R <- signif(cor$estimate,digits = 2)
    
    cor_df[i,] <- c(R,pval)
    
  }
  
  rownames(cor_df) <- key_df$key
  colnames(cor_df) <- c("R","pval")
  
  return(cor_df)
  
}

getBarList <- function(Aclu,Bclu,barcluTable)
{
  AllCluster <- unique(barcluTable$Cluster)
  
  AcluBar <- as.vector(barcluTable$Barcode[which(barcluTable$Cluster == Aclu)])
  BcluBar <- as.vector(barcluTable$Barcode[which(barcluTable$Cluster == Bclu)])
  
  result <- list(AcluBar,BcluBar)
  return(result)
}

#main
GCMat <- readRDS("../data/pat_data.Rdata")

BarCluFile <- "../data/pat.barcodetype.txt"
BarCluTable <- read.table(BarCluFile,sep = "\t",header = TRUE,stringsAsFactors = FALSE)

LigClu <- "B cells"
RecClu <- "Secretory"
cells <- getBarList(RecClu,LigClu,BarCluTable)[[1]]

#get RecTF and TFTar
edge <- read.csv("../net/vis.edge.csv")
Rec <- unique(edge[edge$Source_type == "Rec",1])
TF <- unique(edge[edge$Source_type == "TF",1])
Tar <- unique(edge[edge$Target_type == "Tar",2])

df_RecTF <- unlist(lapply(Rec, function(rec){paste(rec,TF,sep = "_")}))
df_TFTar <- unlist(lapply(TF, function(tf){paste(tf,Tar,sep = "_")}))
df_RecTar <- unlist(lapply(Rec, function(rec){paste(rec,Tar,sep = "_")}))

#get cor of RecTF and TFTar
cor_RecTF <- getCorR(GCMat,df_RecTF,cells)
cor_RecTF <- as.data.frame(cor_RecTF)
cor_RecTF$key <- rownames(cor_RecTF)
cor_RecTF$Rec <- unlist(lapply(cor_RecTF$key, function(x){strsplit(x,"_")[[1]][1]}))
cor_RecTF$TF <- unlist(lapply(cor_RecTF$key, function(x){strsplit(x,"_")[[1]][2]}))
write.table(cor_RecTF,"../cor/cor_RecTF.txt",quote = F,sep = "\t")

cor_TFTar <- getCorR(GCMat,df_TFTar,cells)
cor_TFTar <- as.data.frame(cor_TFTar)
cor_TFTar$key <- rownames(cor_TFTar)
cor_TFTar$TF <- unlist(lapply(cor_TFTar$key, function(x){strsplit(x,"_")[[1]][1]}))
cor_TFTar$Tar <- unlist(lapply(cor_TFTar$key, function(x){strsplit(x,"_")[[1]][2]}))
write.table(cor_TFTar,"../cor/cor_TFTar",quote = F,sep = "\t")

cor_RecTar <- getCorR(GCMat,df_RecTar,cells)
cor_RecTar <- as.data.frame(cor_RecTar)
cor_RecTar$key <- rownames(cor_RecTar)
cor_RecTar$Rec <- unlist(lapply(cor_RecTar$key, function(x){strsplit(x,"_")[[1]][1]}))
cor_RecTar$Tar <- unlist(lapply(cor_RecTar$key, function(x){strsplit(x,"_")[[1]][2]}))
write.table(cor_RecTar,"../cor/cor_RecTar",quote = F,sep = "\t")

#################################################################
#plot
library(ggplot2)
library(ggpubr)
library("RColorBrewer")

display.brewer.pal(n = 8, name = 'OrRd')
mycolor  <- brewer.pal(n = 8, name = 'OrRd')

df_rectf <- cor_RecTF
df_rectf$pval <- -log(df_rectf$pval)
p_RecTF <- ggplot(df_rectf, aes(x = TF,y = Rec)) + 
  geom_point(aes(size=R,color=pval)) + 
  scale_colour_gradient(low=mycolor[4],high=mycolor[8],limits=c(0,175)) + 
  scale_size_continuous(limits = c(-0.1,0.4)) +
  theme_bw() + labs(x = "Transcription factor", y = "Receptor", title = "",color="-1og10(pval)") +
  theme(panel.border = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        legend.position = "bottom",
        legend.direction = "horizontal") 
p_legend <- get_legend(p_RecTF)

df_tftar <- cor_TFTar
df_tftar$pval <- -log(df_tftar$pval)
p_TFTar <- ggplot(df_tftar, aes(x = TF,y = Tar)) + 
  geom_point(aes(size=R,color=pval)) + 
  scale_colour_gradient(low=mycolor[4],high=mycolor[8],limits=c(0,175)) + 
  scale_size_continuous(limits = c(-0.1,0.4)) +
  theme_bw() + labs(x = "", y = "", title = "",color="-1og10(pval)") +  
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line.y = element_line(color = "black")) 

df_rectar <- cor_RecTar
df_rectar$pval <- -log(df_rectar$pval)
p_RecTar <- ggplot(df_rectar, aes(x = Tar,y = Rec)) + 
  geom_point(aes(size=R,color=pval)) + 
  scale_colour_gradient(low=mycolor[4],high=mycolor[8],limits=c(0,175)) + 
  scale_size_continuous(limits = c(-0.1,0.4)) +
  theme_bw() + labs(x = "", y = "", title = "",color="-1og10(pval)") +  
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(color = "black")) 

p_merge <- ggarrange(p_TFTar, NULL, p_RecTF, p_RecTar, 
                ncol = 2, nrow = 2,  align = "hv", legend.grob = p_legend,
                widths = c(4, 1), heights = c(3, 10))
ggsave(p_merge,filename = "../figure/cor_RecTFTar.pdf",width = 8, height = 10)