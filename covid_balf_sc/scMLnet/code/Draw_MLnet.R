DrawMLnet <- function(MLnetList,LigClu,RecClu,workdir,PyHome,plotMLnet = F){
  
  #creat output workdir
  if(!file.exists("output")){ dir.create("output") }
  i <- paste(strsplit(LigClu,split = "\\W")[[1]][1],strsplit(RecClu,split = "\\W")[[1]][1],sep = "_")
  if(is.null(workdir)){
    wd <- paste("./output",i,sep = "/")
    if(!file.exists(wd)){dir.create(wd)}
  }else{
    wd <- paste("./output",workdir,sep = "/")
    if(!file.exists(wd)){dir.create(wd)}
    wd <- paste("./output",workdir,i,sep = "/")
    if(!file.exists(wd)){dir.create(wd)}
  }
  
  #output
  LigRecNet <- MLnetList[[1]]
  RecTFNet <- MLnetList[[2]]
  TFTarNet <- MLnetList[[3]]
  
  #output result
  cat("Save Results\n")
  NetLigRecFile <- paste(wd,"LigRec.net.txt",sep = "/")
  NetRecTFFile <- paste(wd,"RecTF.net.txt",sep = "/")
  NetTFTarFile <- paste(wd,"TFTarGene.net.txt",sep = "/")
  writeLines(LigRecNet,con=NetLigRecFile,sep="\n")
  writeLines(RecTFNet,con=NetRecTFFile,sep="\n")
  writeLines(TFTarNet,con=NetTFTarFile,sep="\n")
  
  #draw MLnet
  if(plotMLnet){
    
    cat("Draw MLnet!\n")
    
    #check python home
    if(is.null(PyHome)){stop("Can't find python!")}
    
    #figure path
    netPic <- paste(wd,"mulnet.pdf",sep = "/")
    
    #draw 
    cmd <- paste(PyHome,"./code/DrawNetNew.py",NetLigRecFile,NetRecTFFile,NetTFTarFile,netPic,sep=" ")
    system(cmd)
    
  }
  
  cat("Finish!\n")
  
}

