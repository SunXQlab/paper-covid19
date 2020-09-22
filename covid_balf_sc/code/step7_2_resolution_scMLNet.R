setwd("./scMLnet")
options(stringsAsFactors = F)

library(Matrix)
library(Seurat)

##Input
GCMat <- readRDS("../data/pat_data.Rdata")
#GCMat<- as(GCMat,"dgCMatrix")

##Source code
source("./code/Run_scMLnet.R")
source("./code/Draw_MLnet.R")


for (res in c("04","08","12","16","2")) {
  
  cat(paste("start running res",res,"!\n",sep = ""))
  
  BarCluFile <- paste("../resolution/res",res,"/res",res,"_barcodetype.txt",sep = "")
  BarCluTable <- read.table(BarCluFile,sep = "\t",header = TRUE,stringsAsFactors = FALSE)
  
  ##get LigClu
  LigClus <- unique(BarCluTable$Cluster)
  LigClus <- LigClus[-grep("Secretory",LigClus)]
  
  ##step1 creat MLnet
  netList <- list()
  for(ligclu in LigClus){
    
    tryCatch({
      
      #ligand cell and receptor cell
      LigClu <- ligclu
      RecClu <- "Secretory cells"
      
      name <- paste(strsplit(LigClu,split = "\\W")[[1]][1],strsplit(RecClu,split = "\\W")[[1]][1],sep = "_")
      print(name)
      
      #main
      netList[[name]] <- CreateMulNet(GCMat,BarCluFile,RecClu,LigClu,res)
      
    }, error=function(e){cat("#--- running error!---#\n\n\n")})
    
  }
  
  #step2 save output and plot MLnet
  #output path
  workdir <- paste("res",res,sep = "")
  #python home
  PyHome <- "D:/Miniconda3/envs/R36/python.exe"
  for (name in names(netList)) {
    
    #MLnet
    MLnetList <- netList[[name]]
    print(paste0(name,":"))
    
    #ligand cell and receptor cell
    LigClu <- strsplit(name,"_")[[1]][1]
    RecClu <- strsplit(name,"_")[[1]][2]
    
    #only save output
    DrawMLnet(MLnetList,LigClu,RecClu,workdir,plotMLnet = F)
    
  }
  
  cat(paste("finish running res",res,"!\n\n\n",sep = ""))
  
}


