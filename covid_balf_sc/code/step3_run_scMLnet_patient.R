setwd("./scMLnet")
options(stringsAsFactors = F)

library(Matrix)
library(Seurat)
 
##Input
GCMat <- readRDS("../data/pat_data.Rdata")
#GCMat<- as(GCMat,"dgCMatrix")

BarCluFile <- "../data/pat.barcodetype.txt"
BarCluTable <- read.table(BarCluFile,sep = "\t",header = TRUE,stringsAsFactors = FALSE)

##get LigClu
LigClus <- unique(BarCluTable$Cluster)
#skip Secretory(RecClu) and Doublets(low quality)
LigClus <- LigClus[-grep("Secretory|Doublets",LigClus)]

##Parameter(Default)
# pval <- 0.05
# logfc <- 0.15
# LigRecLib <- "./database/LigRec.txt"
# TFTarLib <- "./database/TFTargetGene.txt"
# RecTFLib <- "./database/RecTF.txt"

##Source code
source("./code/Run_scMLnet.R")
source("./code/Draw_MLnet.R")

###################################################################

##step1 creat MLnet
netList <- list()
for(ligclu in LigClus){
  
  tryCatch({
    
    #ligand cell and receptor cell
    LigClu <- ligclu
    RecClu <- "Secretory"
    
    name <- paste(strsplit(LigClu,split = "\\W")[[1]][1],RecClu,sep = "_")
    
    #main
    netList[[name]] <- CreateMulNet(GCMat,BarCluFile,RecClu,LigClu)
    
  }, error=function(e){cat("#--- running error!---#\n")})
  
}

#step2 save output and plot MLnet
#output path
workdir <- "nCoV_patient"
for (name in names(netList)) {
  
  
  #MLnet
  MLnetList <- netList[[name]]
  print(paste0(name,":"))
  
  #ligand cell and receptor cell
  LigClu <- strsplit(name,"_")[[1]][1]
  RecClu <- strsplit(name,"_")[[1]][2]
  
  #main
  #only save output
  DrawMLnet(MLnetList,LigClu,RecClu,workdir,plotMLnet = F)
  
}

