#install.packages("VennDiagram")
library(VennDiagram)
library(ggpubr)

setwd("scMLnet/")
options(stringsAsFactors = F)

###############

#node

##input

wds <- paste0("res",c("04","08","12","16","2"),"_ACE2")

Ligand=list()
Receptor=list()
TF=list()

for(wd in wds){
  
  wd_in <- paste0("./output/",wd,"/")
  folder.list <- list.files(wd_in)
  
  lig=c()
  rec=c()
  tf=c()
  
  for(folder in folder.list){
    
    print(paste(wd,folder,sep = ":"))
    print(folder)
    workdir <- paste0(wd_in,folder)
    file.list <- list.files(workdir)
    
    LigRecTable <- read.table(paste(workdir,"LigRec.net.txt",sep = "/"))
    lig_tmp <- unlist(lapply(LigRecTable$V1,function(x){strsplit(x,"_")[[1]][1]}))
    rec_tmp <- unlist(lapply(LigRecTable$V1,function(x){strsplit(x,"_")[[1]][2]}))
    
    RecTFTable <- read.table(paste(workdir,"RecTF.net.txt",sep = "/"))
    tf_tmp <- unlist(lapply(RecTFTable$V1,function(x){strsplit(x,"_")[[1]][2]}))
    
    lig=c(lig,lig_tmp)
    rec=c(rec,rec_tmp)
    tf=c(tf,tf_tmp)
    
  }
  
  name <- unlist(lapply(wd, function(x){strsplit(x,split = "_")[[1]][1]}))
  Ligand[[name]]=unique(lig)
  Receptor[[name]]=unique(rec)
  TF[[name]]=unique(tf)
  
}
names(Ligand) <- paste0("res=",c("0.4","0.8","1.2","1.6","2"))
names(Receptor) <- paste0("res=",c("0.4","0.8","1.2","1.6","2"))
names(TF) <- paste0("res=",c("0.4","0.8","1.2","1.6","2"))

ls_mer <- list(Ligand=Ligand,
               Receptor=Receptor,
               TF=TF)

##plot
###venn
for (i in c("Ligand","Receptor","TF")) {
  
  data <- ls_mer[[i]]
  
  p_venn <- venn.diagram(data,filename=NULL, cex = 1.5,
                         cat.cex = 1.5,
                         col = "black",fill = c("cornflowerblue", "green", "yellow", "darkorchid1","LightPink"),
                         alpha = 0.4,cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4","PaleVioletRed"))
  
  dpi=500
  pdf(file=paste0("../resolution/venn_overlap_",i,".pdf"), width = 6, height = 6)
  print(grid.draw(p_venn))
  dev.off()
  
}

###############

#links

##input

wds <- paste0("res",c("04","08","12","16","2"),"_ACE2")

LigRec=list()
RecTF=list()

for(wd in wds){
  
  wd_in <- paste0("./output/",wd,"/")
  folder.list <- list.files(wd_in)
  
  ligrec=c()
  rectf=c()
  
  for(folder in folder.list){
    
    print(paste(wd,folder,sep = ":"))
    print(folder)
    workdir <- paste0(wd_in,folder)
    file.list <- list.files(workdir)
    
    LigRecTable <- read.table(paste(workdir,"LigRec.net.txt",sep = "/"))
    ligrec_tmp <- LigRecTable$V1
    
    RecTFTable <- read.table(paste(workdir,"RecTF.net.txt",sep = "/"))
    rectf_tmp <- RecTFTable$V1
    
    ligrec=c(ligrec,ligrec_tmp)
    rectf=c(rectf,rectf_tmp)
    
  }
  
  name <- unlist(lapply(wd, function(x){strsplit(x,split = "_")[[1]][1]}))
  LigRec[[name]]=unique(ligrec)
  RecTF[[name]]=unique(rectf)
  
}
names(LigRec) <- paste0("res=",c("0.4","0.8","1.2","1.6","2"))
names(RecTF) <- paste0("res=",c("0.4","0.8","1.2","1.6","2"))

ls_mer <- list(LigRec=LigRec,
               RecTF=RecTF)

##plot
###venn
for (i in c("LigRec","RecTF")) {
  
  data <- ls_mer[[i]]
  
  p_venn <- venn.diagram(data,filename=NULL, cat.cex = 1.5, cex = 1.5,
                         col = "black",fill = c("cornflowerblue", "green", "yellow", "darkorchid1","LightPink"),
                         alpha = 0.4,cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4","PaleVioletRed"));
  
  pdf(file=paste0("../resolution/venn_overlap_pair_",i,".pdf"), width = 6, height = 6)
  print(grid.draw(p_venn))
  dev.off()
  
}

###############

#count

##input
wds <- paste0("res",c("04","08","12","16","2"),"_ACE2")

Ligand=list()
Receptor=list()
TF=list()

for(wd in wds){
  
  wd_in <- paste0("./output/",wd,"/")
  folder.list <- list.files(wd_in)
  
  lig=c()
  rec=c()
  tf=c()
  
  for(folder in folder.list){
    
    print(paste(wd,folder,sep = ":"))
    print(folder)
    workdir <- paste0(wd_in,folder)
    file.list <- list.files(workdir)
    
    LigRecTable <- read.table(paste(workdir,"LigRec.net.txt",sep = "/"))
    lig_tmp <- unlist(lapply(LigRecTable$V1,function(x){strsplit(x,"_")[[1]][1]}))
    rec_tmp <- unlist(lapply(LigRecTable$V1,function(x){strsplit(x,"_")[[1]][2]}))
    
    RecTFTable <- read.table(paste(workdir,"RecTF.net.txt",sep = "/"))
    tf_tmp <- unlist(lapply(RecTFTable$V1,function(x){strsplit(x,"_")[[1]][2]}))
    
    lig=c(lig,lig_tmp)
    rec=c(rec,rec_tmp)
    tf=c(tf,tf_tmp)
    
  }
  
  name <- unlist(lapply(wd, function(x){strsplit(x,split = "_")[[1]][1]}))
  Ligand[[name]]=unique(lig)
  Receptor[[name]]=unique(rec)
  TF[[name]]=unique(tf)
  
}
names(Ligand) <- paste0("res=",c("0.4","0.8","1.2","1.6","2"))
names(Receptor) <- paste0("res=",c("0.4","0.8","1.2","1.6","2"))
names(TF) <- paste0("res=",c("0.4","0.8","1.2","1.6","2"))

count_node <- data.frame(resolution=c("0.4","0.8","1.2","1.6","2"),
                         Ligand=lengths(Ligand),
                         Receptor=lengths(Receptor),
                         TF=lengths(TF))
count_df <- melt(count_node,id.vars = "resolution")

##plot
p1 <- ggplot(count_df,aes(x=resolution,y=value)) + 
  geom_line(aes(col=variable,group=variable),size=1) +
  theme_bw() + labs(x = "resolution")+   
  coord_cartesian(ylim = c(0,5)) +   
  scale_y_continuous(breaks = c(0,5,5)) +  
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 24),
        axis.text = element_text(size = 22),
        axis.line = element_line(color = "black"),
        legend.position = "bottom",
        legend.direction = "horizontal") 

p2 <- ggplot(count_df,aes(x=resolution,y=value)) + 
  geom_line(aes(col=variable,group=variable),size=1) +
  theme_bw() + labs(x = "", y = "count")+
  coord_cartesian(ylim = c(45,65)) +  
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 22),
        axis.title.y = element_text(hjust = 0.4,size = 24),
        legend.title = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.line.y = element_line(color = "black")) 

p_node <- ggarrange(p2,p1,heights=c(4/5,1/4),
                    ncol = 1, nrow = 2, align = "v",
                    common.legend = F) 

pdf(file=paste0("../resolution/count_node.pdf"), width = 8, height = 9)
print(p_node)
dev.off()