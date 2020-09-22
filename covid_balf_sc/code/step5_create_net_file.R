rm(list = ls())
gc()

options(stringsAsFactors = F)
library(scater)

#Function
##calculate aver expr
calaverexpr <- function(data.use,cells){
  
  cpm <- calculateCPM(data.use)
  cpm <- log1p(cpm)
  averexpr <- rowMeans(cpm[, cells])
  
  
  return(averexpr)
  
}

##get Baecode
getBarList <- function(Aclu,data.use,barCluTable){
  
  AllCluster <- unique(barCluTable$Cluster)
  AcluBar <- as.vector(barCluTable$Barcode[which(barCluTable$Cluster == Aclu)])
  return(AcluBar)
  
}

#input path
wd_in <- paste0("./output/nCoV_patient/")
folder.list <- list.files(wd_in)

#output path
wd_out <- paste0("../net")
if(!file.exists(wd_out)){dir.create(wd_out)}

#############################################

#average expr
GCMat <- readRDS("../data/pat_data.Rdata")
BarCluTable <- read.table("../data/pat.barcodetype.txt",sep = "\t",header = TRUE,stringsAsFactors = FALSE)
df_averexpr <- data.frame()
for(clu in unique(BarCluTable$Cluster)){
  
  bar <- getBarList(clu,GCMat,BarCluTable)
  averexpr <- calaverexpr(GCMat,bar)
  averexpr_tmp <- data.frame(gene = rownames(GCMat),
                             averexpr = averexpr,
                             Cluster = rep(clu,nrow(GCMat)))
  df_averexpr <- rbind(df_averexpr,averexpr_tmp)
  
}
write.table(df_averexpr,file = paste(wd_out,"averexpr.txt",sep = "/"),quote = F,sep = "\t",col.names = T,row.names = F)

## get ACE2 Lig-Rec MLnet 
Clusters <- unique(BarCluTable$Cluster)
LigRec_sel <- data.frame()
for(folder in folder.list){
  
  print(folder)
  workdir <- paste0(wd_in,folder)
  file.list <- list.files(workdir)
  
  LigRecTable <- read.table(paste(workdir,"LigRec.net.txt",sep = "/"))
  LigRecTable$Lig <- unlist(lapply(LigRecTable$V1,function(x){strsplit(x,"_")[[1]][1]}))
  LigRecTable$Rec <- unlist(lapply(LigRecTable$V1,function(x){strsplit(x,"_")[[1]][2]}))
  LigRecTable$Type <- rep(folder,nrow(LigRecTable))
  
  LigClu <- Clusters[grep(strsplit(folder,"_")[[1]][1],Clusters)]
  RecClu <- Clusters[grep(strsplit(folder,"_")[[1]][2],Clusters)]
  
  LigRecTable$Score <- unlist(apply(LigRecTable,1,function(LigRec_temp){
    
    lig <- LigRec_temp["Lig"]
    lig_s <- df_averexpr$averexpr[df_averexpr$gene == lig & df_averexpr$Cluster == LigClu]
    
    rec <- LigRec_temp["Rec"]
    rec_s <- df_averexpr$averexpr[df_averexpr$gene == rec & df_averexpr$Cluster == RecClu]
    
    lig_s*rec_s
    
  }))
  
  LigRec_sel <- rbind(LigRec_sel,LigRecTable)
  
}

LigRec_last <- LigRec_sel[order(LigRec_sel$V1,LigRec_sel$Score,decreasing = T),]
LigRec_last <- LigRec_last[!duplicated(LigRec_last$V1),]
LigRec_last <- LigRec_last[LigRec_ord$Score > 1.5,]
colnames(LigRec_last) <- c("Key","Lig","Rec","Type","Score")

write.csv(LigRec_last,paste(wd_out,"vis.LigRec.net.csv",sep = "/"),quote = F,row.names = F)


##################################

RecTF_sel <- data.frame()
for(folder in folder.list){
  
  print(folder)
  workdir <- paste0(wd_in,folder)
  file.list <- list.files(workdir)
  
  RecTFTable <- read.table(paste(workdir,"RecTF.net.txt",sep = "/"))
  RecTFTable$Rec <- unlist(lapply(RecTFTable$V1,function(x){strsplit(x,"_")[[1]][1]}))
  RecTFTable$TF <- unlist(lapply(RecTFTable$V1,function(x){strsplit(x,"_")[[1]][2]}))
  RecTFTable$Type <- rep(folder,nrow(RecTFTable))
  
  Rec_sel <- unique(LigRec_last$Rec[LigRec_last$Type == folder])
  RecTFTable <- RecTFTable[RecTFTable$Rec %in% Rec_sel,]
  
  RecTF_sel <- rbind(RecTF_sel,RecTFTable)
  
}

RecTF_last <- RecTF_sel
colnames(RecTF_last) <- c("Key","Rec","TF","Type")
write.csv(RecTF_last,paste(wd_out,"vis.RecTF.net.csv",sep = "/"),quote = F,row.names = F)

##################################

TFTar_sel <- data.frame()
for(folder in folder.list){
  
  print(folder)
  workdir <- paste0(wd_in,folder)
  file.list <- list.files(workdir)
  
  TFTarTable <- read.table(paste(workdir,"TFTarGene.net.txt",sep = "/"))
  TFTarTable$TF <- unlist(lapply(TFTarTable$V1,function(x){strsplit(x,"_")[[1]][1]}))
  TFTarTable$Tar <- unlist(lapply(TFTarTable$V1,function(x){strsplit(x,"_")[[1]][2]}))
  TFTarTable$Type <- rep(folder,nrow(TFTarTable))
  
  TF_sel <- unique(RecTF_last$TF[RecTF_last$Type == folder])
  TFTarTable <- TFTarTable[TFTarTable$TF %in% TF_sel,]
  
  TFTar_sel <- rbind(TFTar_sel,TFTarTable)
  
}

TFTar_last <- TFTar_sel
colnames(TFTar_last) <- c("Key","TF","Tar","Type")
write.csv(TFTar_last,paste(wd_out,"vis.TFTar.net.csv",sep = "/"),quote = F,row.names = F)

##################################
library(dplyr)

##create node
Lig <- LigRec_last$Lig[]
Rec_uni <- unique(LigRec_last$Rec)
TF_uni <- unique(TFTar_last$TF)
Tar_uni <- unique(TFTar_last$Tar)

node <- data.frame(Node_name=c(Lig,Rec_uni,TF_uni,Tar_uni))
node$node_type <- c(rep("Lig",length(Lig)),rep("Rec",length(Rec_uni)),rep("TF",length(TF_uni)),rep("Tar",length(Tar_uni)))
node$Type <- c(LigRec_last$Type,rep("Secretory",times=nrow(node)-nrow(LigRec_last)))

node <- node[!duplicated(node),]
node <- node[order(node$node_type,node$Node_name,node$Type),]
node$ID <- paste("node",1:nrow(node),sep = ".")
write.csv(node,paste(wd_out,"vis.node.csv",sep = "/"),quote = F,row.names = F)

##create edge
edge_LigRec <- LigRec_last[,c(2,3,1,4)]
edge_LigRec$Source_type <- rep("Lig",nrow(edge_LigRec))
edge_LigRec$Target_type <- rep("Rec",nrow(edge_LigRec))
colnames(edge_LigRec) <- c("Source","Target","Key","Type","Source_type","Target_type")

edge_RecTF <- RecTF_last[,c(2,3,1,4)]
edge_RecTF$Source_type <- rep("Rec",nrow(edge_RecTF))
edge_RecTF$Target_type <- rep("TF",nrow(edge_RecTF))
colnames(edge_RecTF) <- c("Source","Target","Key","Type","Source_type","Target_type")

edge_TFTar <- TFTar_last[,c(2,3,1,4)]
edge_TFTar$Source_type <- rep("TF",nrow(edge_TFTar))
edge_TFTar$Target_type <- rep("Tar",nrow(edge_TFTar))
colnames(edge_TFTar) <- c("Source","Target","Key","Type","Source_type","Target_type")

edge <- rbind(edge_LigRec,edge_RecTF,edge_TFTar)
edge$Source_ID <- unlist(apply(edge,1,function(x){
  
  if(x[5] == "Lig"){
    ID=node$ID[intersect(which(node$Node_name == x[1]),which(node$Type == x[4]))]
  }else{
    ID=node$ID[which(node$Node_name == x[1])]
  }
  ID
  
}))
edge$Target_ID <- unlist(apply(edge,1,function(x){ID=node$ID[which(node$Node_name == x[2])]}))
write.csv(edge,paste(wd_out,"vis.edge.csv",sep = "/"),quote = F,row.names = F)
