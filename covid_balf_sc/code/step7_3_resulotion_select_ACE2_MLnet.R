setwd("./scMLnet")
options(stringsAsFactors = F)

source("./code/Draw_MLnet.R")
for(res in c("04","08","12","16","2")){
  
  workdir_in <- paste0("./output/res",res,"/")
  folder.list <- list.files(workdir_in)
  
  for(folder in folder.list){
    
    print(folder)
    wd <- paste0(workdir_in,folder)
    file.list <- list.files(wd)
    
    #read MLnet
    LigRecfile <- file.list[grep("LigRec.net.txt",file.list)]
    TFTarfile <- file.list[grep("TFTarGene.net.txt",file.list)]
    RecTFfile <- file.list[grep("RecTF.net.txt",file.list)]
    
    LigRec <- read.table(paste(wd,LigRecfile,sep = "/"))
    LigRec$Lig <- unlist(lapply(LigRec$V1,function(x){strsplit(x,"_")[[1]][1]}))
    LigRec$Rec <- unlist(lapply(LigRec$V1,function(x){strsplit(x,"_")[[1]][2]}))
    
    RecTF <- read.table(paste(wd,RecTFfile,sep = "/"))
    RecTF$Rec <- unlist(lapply(RecTF$V1,function(x){strsplit(x,"_")[[1]][1]}))
    RecTF$TF <- unlist(lapply(RecTF$V1,function(x){strsplit(x,"_")[[1]][2]}))
    
    TFTar <- read.table(paste(wd,TFTarfile,sep = "/"))
    TFTar$TF <- unlist(lapply(TFTar$V1,function(x){strsplit(x,"_")[[1]][1]}))
    TFTar$Tar <- unlist(lapply(TFTar$V1,function(x){strsplit(x,"_")[[1]][2]}))
    
    ##Creat ACE2 MLnet
    TFTar <- TFTar[grep("ACE2",TFTar$Tar),]
    RecTF <- RecTF[RecTF$TF %in% unique(TFTar$TF),]
    LigRec <- LigRec[LigRec$Rec %in% unique(RecTF$Rec),]
    
    if(nrow(TFTar) & nrow(RecTF) & nrow(LigRec)){
      
      #Draw and Save ACE2 MLnet
      MLnetList <- list("LigRec" = LigRec$V1,
                        "RecTF" = RecTF$V1,
                        "TFTar" = TFTar$V1)
      workdir_out <- paste("res",res,"_ACE2",sep = "")
      LigClu <- strsplit(folder,"_")[[1]][1]
      RecClu <- strsplit(folder,"_")[[1]][2]
      PyHome <- "D:/Miniconda3/envs/R36/python.exe"
      DrawMLnet(MLnetList,LigClu,RecClu,workdir_out,PyHome,plotMLnet = T)
      
    }
    
  }
  
}


