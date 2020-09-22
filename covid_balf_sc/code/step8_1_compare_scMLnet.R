setwd("./compare")
library(igraph)

options(stringsAsFactors = F)

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

geneset =  read.table("./DEG/geneset.txt",header = F,sep = "\t")
geneset <- unlist(geneset)

workdir_in <- paste0("../scMLnet/output/nCoV_patient/")
folder.list <- list.files(workdir_in)

df_RecTF <- c()
df_TFTar <- c()
df_RecTar <- c()
df_LigTar <- c()
for(folder in folder.list){
  
  print(folder)
  wd <- paste0(workdir_in,folder)
  file.list <- list.files(wd)
  
  #read MLnet
  TFTarfile <- file.list[grep("TFTarGene.net.txt",file.list)]
  RecTFfile <- file.list[grep("RecTF.net.txt",file.list)]
  LigRecfile <- file.list[grep("LigRec.net.txt",file.list)]
  
  TFTar <- read.table(paste(wd,TFTarfile,sep = "/"))
  TFTar$from <- unlist(lapply(TFTar$V1,function(x){strsplit(x,"_")[[1]][1]}))
  TFTar$to <- unlist(lapply(TFTar$V1,function(x){strsplit(x,"_")[[1]][2]}))
  
  TFTar <- TFTar[TFTar$to %in% geneset,]
  Tar <- unique(TFTar$to)
  TF <- unique(TFTar$from)
  
  RecTF <- read.table(paste(wd,RecTFfile,sep = "/"))
  RecTF$from <- unlist(lapply(RecTF$V1,function(x){strsplit(x,"_")[[1]][1]}))
  RecTF$to <- unlist(lapply(RecTF$V1,function(x){strsplit(x,"_")[[1]][2]}))
  
  RecTF <- RecTF[RecTF$to %in% TF,]
  Rec <- unique(RecTF$from)
  
  LigRec <- read.table(paste(wd,LigRecfile,sep = "/"))
  LigRec$from <- unlist(lapply(LigRec$V1,function(x){strsplit(x,"_")[[1]][1]}))
  LigRec$to <- unlist(lapply(LigRec$V1,function(x){strsplit(x,"_")[[1]][2]}))
  
  LigRec <- LigRec[LigRec$to %in% Rec,]
  Lig <- unique(LigRec$from)
  
  df_RecTF <- c(df_RecTF,RecTF$V1)
  df_TFTar <- c(df_TFTar,TFTar$V1)
  
  net <- rbind(LigRec[,2:3],RecTF[,2:3],TFTar[,2:3])
  g <- graph_from_data_frame(net, directed=TRUE)
  d  <- shortest.paths(g,mode="out",weights = NA)
  
  is.RecTar.link <- c()
  for (i in Rec) {
    for (j in Tar) {
      is.RecTar.link <- c(is.RecTar.link,!is.infinite(d[i,j]))
    }
  }
  names(is.RecTar.link) <- paste(rep(Rec,each=length(Tar)),Tar,sep = "_")
  
  df_RecTar <- c(df_RecTar,names(is.RecTar.link)[is.RecTar.link])
  
  is.LigTar.link <- c()
  for (i in Lig) {
    for (j in Tar) {
      is.LigTar.link <- c(is.LigTar.link,!is.infinite(d[i,j]))
    }
  }
  names(is.LigTar.link) <- paste(rep(Lig,each=length(Tar)),Tar,sep = "_")
  
  df_LigTar <- c(df_LigTar,names(is.LigTar.link)[is.LigTar.link])
  
}

df_TFTar <- unique(df_TFTar)
df_RecTF <- unique(df_RecTF)
df_RecTar <- unique(df_RecTar)

df_LigTar <- unique(df_LigTar)
write.table(df_LigTar,"./scMLnet/df_LigTar.txt",quote = F,sep = "\t",row.names = F,col.names = F)

#get cor of RecTF and TFTar
cor_RecTF <- getCorR(GCMat,df_RecTF,cells)
write.csv(cor_RecTF,file = "./scMLnet/cor_RecTF.csv")

cor_TFTar <- getCorR(GCMat,df_TFTar,cells)
write.csv(cor_TFTar,file = "./scMLnet/cor_TFTar.csv")

cor_RecTar <- getCorR(GCMat,df_RecTar,cells)
write.csv(cor_RecTar,file = "./scMLnet/cor_RecTar.csv")
