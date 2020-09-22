options(stringsAsFactors = F)

#Function
Fisher_test <- function(gene_list,gene_up,TotalGene){
  
  a <- length(intersect(gene_list,gene_up))
  b <- length(gene_list)-a
  c <- length(gene_up)-a
  d <- length(TotalGene)-a-b-c
  
  matrix <- matrix(c(a,c,b,d),nrow=2)
  res <- fisher.test(matrix,alternative="greater")
  
  return(res$p.value)
}

#Input

##DEG and norm data
deg=read.table('./output/DE_genes.txt',sep = '\t',header = T)
deg=na.omit(deg)
data=read.table("./data/norm.txt",sep = '\t',header = T)
boxplot(log1p(data))

gene_mean=apply(data,1,mean)
deg_up=deg$gene[deg$log2FoldChange > log2(4) & deg$padj < 0.01]
gene_expressed=as.matrix(rownames(data)[gene_mean > 10])

#database

##read ligand receptor information

LigRec=read.delim("database/LigRec.txt")

LigRec=toupper(as.matrix(LigRec))
LigRec=LigRec[,2:3]

Ligand=unique(LigRec[,1])
Receptor=unique(LigRec[,2])

##read TF Target information
TF_target=read.table("database/TFTargetGene.txt",header = T,sep="\t")
TF_target=as.matrix(TF_target)
TF_target=TF_target[,-3]


##read KEGG pathway information (graphite package)

KEGG_all_edge=read.delim("database/KEGG_all_edge.txt")
KEGG_all_edge=as.matrix(KEGG_all_edge)
KEGG_pathway=unique(KEGG_all_edge[,5])
KEGG_edge_symbol=list()
for (i in 1:length(KEGG_pathway)){
  KEGG_edge_symbol[[i]]=KEGG_all_edge[KEGG_all_edge[,5]%in%KEGG_pathway[i],]
}
names(KEGG_edge_symbol)=KEGG_pathway

##read heterodimer information

KEGG_binding=c("Cytokine-cytokine receptor interaction","ECM-receptor interaction","Neuroactive ligand-receptor interaction")
KEGG_binding=do.call(rbind,KEGG_edge_symbol[KEGG_binding])
KEGG_binding=KEGG_binding[KEGG_binding[,4]=="Binding",]
KEGG_binding=rbind(KEGG_binding[,1:2],KEGG_binding[,2:1])

Manual_binding=read.delim("database/manual_binding.txt",header=F)
Manual_binding=as.matrix(Manual_binding)
Manual_binding=rbind(Manual_binding,Manual_binding[,2:1])

Receptor_binding=rbind(KEGG_binding,Manual_binding)
Receptor_binding=Receptor_binding[!duplicated(Receptor_binding),]

#main

##part1

###get KEGG Ligand-Receptor info

KEGG_LR=merge(KEGG_all_edge,LigRec,by.x=1:2,by.y=1:2,sort=F)
KEGG_LR=as.matrix(KEGG_LR)

expressed_ligand=intersect(gene_expressed,Ligand)
expressed_receptor=intersect(gene_expressed,Receptor)

LR_select=KEGG_LR[KEGG_LR[,1]%in%expressed_ligand & KEGG_LR[,2]%in%expressed_receptor,]

KEGG_LR_select=list()
for (i in 1:length(KEGG_pathway)){
  KEGG_LR_select[[i]]=LR_select[LR_select[,5]%in%KEGG_pathway[i],]
}
names(KEGG_LR_select)=KEGG_pathway

###get KEGG TF info

KEGG_TF=list()
for (i in 1:length(KEGG_edge_symbol)){
  tmp=KEGG_edge_symbol[[i]]
  if(!is.matrix(tmp)){tmp=t(as.matrix(tmp))}
  KEGG_TF[[i]]=intersect(unique(as.vector(as.matrix(tmp[,1:2]))),TF_target[,1])
}
names(KEGG_TF)=names(KEGG_edge_symbol)

###find activated transcription factors

TF_list=unique(TF_target[,1])
TF_fisher=c()
for (i in 1:length(TF_list))
{
  x=TF_target[TF_target[,1]==TF_list[i],2]
  TF_fisher[i]=Fisher_test(gene_list = x, gene_up = deg_up, TotalGene = as.matrix(rownames(data)))
}
names(TF_fisher)=TF_list
TF_diff=TF_list[TF_fisher<0.05]
TF_diff=intersect(gene_expressed[,1],TF_diff)

write.table(data.frame(TF=TF_diff,
                       pval=as.character(TF_fisher[TF_diff])),
            file = "output/fisher_TF.txt",sep = "\t",quote = F)

###select pathways by activated transcription factors

KEGG_TF_select=list()
for (i in 1:length(KEGG_edge_symbol)){
  x=KEGG_TF[[i]]
  KEGG_TF_select[[i]]=intersect(KEGG_TF[[i]],TF_diff)
}
names(KEGG_TF_select)=names(KEGG_edge_symbol)

###find activated pathways

pathway_select=names(KEGG_edge_symbol)[unlist(lapply(KEGG_LR_select,length)>0) &
                                         unlist(lapply(KEGG_TF_select,length)>0)]

##part2

###get LR and TFs in activated pathways

pathway_select_LR_list=KEGG_LR_select[pathway_select]
pathway_select_TF_list=KEGG_TF_select[pathway_select]

###get heterodimer LR in activated pathways

x=pathway_select_LR_list
pathway_heterodimer=list()
for (i in 1:length(x)){
  y=x[[i]]
  if (is.matrix(y)){y=y[!duplicated(y),]}
  if (!is.matrix(y)){y=t(as.matrix(y))}
  pathway_heterodimer[[i]]=matrix(ncol=2)
  if(nrow(y)>0){
    for (j in 1:nrow(y)){
      y1=y[j,1]
      y2=y[j,2]
      y3=unique(as.matrix(merge(y2,Receptor_binding,by.x=1,by.y=1,sort=F)[,2]))
      y4=merge(y3,LigRec,by.x=1,by.y=2,sort=F)
      y4=as.matrix(y4[,2:1])
      y5=as.matrix(merge(y1,y4,by.x=1,by.y=1,sort=F))
      if (length(y5)>0){
        y6=cbind(y2,y5[,2])
        pathway_heterodimer[[i]]=rbind(pathway_heterodimer[[i]],y6)
      }
    }
  }else{
    pathway_heterodimer[[i]]=pathway_heterodimer[[i]]
  }
}
pathway_heterodimer=lapply(pathway_heterodimer,function(x){x=x[-1,]})

pathway_select_LR_list=x
names(pathway_heterodimer)=names(x)

###get sub-networks in activated pathways

library("igraph")
x=pathway_select
y=KEGG_edge_symbol[x]

pathway_branch=list()
for (i in 1:length(y)){
  yy=y[[i]][,1:2]
  b=pathway_select_LR_list[[i]]
  if(!is.matrix(b)){b=t(as.matrix(b))}
  yy2=pathway_heterodimer[[i]]
  if (!is.matrix(yy2)){yy2=t(as.matrix(yy2))}
  
  yyy=yy[yy[,2]%in%Receptor&(!yy[,2]%in%b[,2])&(!yy[,2]%in%as.vector(yy2)),2] 
  yy=yy[!(yy[,2]%in%yyy|yy[,1]%in%yyy),]
  yy=yy[!(yy[,1]%in%Ligand|yy[,2]%in%Ligand),] 
  yy=rbind(yy,b[,1:2])
  
  yy=rbind(yy,yy2) 
  yy=rbind(yy,yy2[,2:1]) 
  yy=yy[!duplicated(yy),] 
  yy=graph.edgelist(yy) 
  a=shortest.paths(yy,mode="out") 
  b1=unique(b[,2])
  d=pathway_select_TF_list[[i]]
  g=a[b1,] 
  if (!is.matrix(g)){
    g=t(as.matrix(g))
    rownames(g)=b1
  }
  g1=apply(g,1,function(x){colnames(a)[x<50]})
  g1=unique(as.vector(unlist(g1)))
  g2=a[g1,] 
  if (!is.matrix(g2)){
    g2=t(as.matrix(g2))
    rownames(g2)=g1
  }
  g3=apply(g2,1,function(x){y=x[match(d,colnames(a))];y=names(y)[y<10000]})
  g4=lapply(g3,length)
  g5=names(g4)[g4>0] 
  pathway_branch[[i]]=g5 
}
names(pathway_branch)=x[1:length(pathway_branch)]

###get edge of sub-networks in activated pathways

x=pathway_branch
y=lapply(x,is.null)
x=x[!unlist(y)]

LR_x=pathway_select_LR_list[names(x)]

pathway_sub_net=list()
for (i in 1:length(x)){
  y=KEGG_edge_symbol[names(x)[i]][[1]][,1:4]
  
  
  yy2=pathway_heterodimer[names(x)][[i]]
  if (!is.matrix(yy2)){yy2=t(as.matrix(yy2))}
  tmp=LR_x[[i]]
  check=is.matrix(tmp)
  if (!check){tmp=t(as.matrix(tmp))}
  yy=y[y[,2]%in%Receptor&(!y[,2]%in%tmp[,2])&(!y[,2]%in%as.vector(yy2)),2]
  yy=y[!(y[,2]%in%yy|y[,1]%in%yy),]
  yy=yy[!(yy[,1]%in%Ligand|yy[,2]%in%Ligand),]
  yy=yy[!(yy[,2]%in%Receptor),]
  yy=yy[!duplicated(yy),]
  y=yy
  
  y1=y[y[,1]%in%x[[i]]&y[,2]%in%x[[i]],]
  if (check) {y2=cbind(tmp[,1:2],"directed","LR")}
  if (!check){y2=t(as.matrix(c(tmp[,1:2],"directed","LR")))}
  y3=pathway_heterodimer[names(x)][[i]]
  if (!is.matrix(y3)){y3=t(as.matrix(y3))}
  y3=cbind(y3,"directed","heterodimer")
  y2=y2[y2[,2]%in%y1[,1]|y2[,2]%in%y3[,1],]
  if (!is.matrix(y2)){y2=t(as.matrix(y2))}
  colnames(y2)=c("From","To","Direction","Interaction_type")
  pathway_sub_net[[i]]=rbind(y2,y3,y1)
}
names(pathway_sub_net)=names(x)

###get node of sub-networks in activated pathways

pathway_sub_node=list()
x=names(pathway_sub_net)
x1=pathway_select_LR_list[x]
x2=pathway_select_TF_list[x]
for (i in 1:length(x)){
  y=pathway_sub_net[[i]]
  y=unique(as.vector(y[,1:2]))
  y=cbind(y,"links")
  yy=x1[[i]]
  check=is.matrix(yy)
  if (!check){yy=t(as.matrix(yy))}
  y[y[,1]%in%x2[[i]],2]="TF"
  y[y[,1]%in%yy[,2],2]="Receptor"
  y[y[,1]%in%yy[,1],2]="Ligand"
  colnames(y)=c("gene","gene_type")
  pathway_sub_node[[i]]=y
}
names(pathway_sub_node)=names(pathway_sub_net)

###find activated pathway branch

x=pathway_branch
y=lapply(x,is.null)
x=x[!unlist(y)]

TF_x=pathway_select_TF_list[names(x)]
TF_x=unique(unlist(TF_x))

Pathway_fisher=c()
for (i in 1:length(x))
{
  Pathway_fisher[i] = Fisher_test(gene_list=x[[i]], gene_up = unique(c(as.matrix(deg_up),TF_x)), TotalGene = as.matrix(rownames(data)))
}
names(Pathway_fisher)=names(x)
Pathway_fisher=as.matrix(Pathway_fisher)
Pathway_fisher=Pathway_fisher[order(Pathway_fisher[,1]),]

write.table(data.frame(pathway=names(Pathway_fisher),
                       pval=as.character(Pathway_fisher)),
            file = "output/fisher_pathway.txt",
            sep = "\t",quote = F)

#Output

##remove pathway with p<0.05

filter_test=names(Pathway_fisher)[Pathway_fisher<0.05]
pathway_branch_final_network=pathway_sub_net[filter_test]
pathway_branch_final_node_attributes=pathway_sub_node[filter_test]

pathway_combine=do.call(rbind,pathway_branch_final_network)
pathway_combine=pathway_combine[!duplicated(pathway_combine[,1:2]),]

pathway_combine_node=do.call(rbind,pathway_branch_final_node_attributes)
tmp1=which(pathway_combine_node[,2]=="Ligand")
tmp2=which(pathway_combine_node[,2]=="Receptor")
tmp3=which(pathway_combine_node[,2]=="links")
tmp4=which(pathway_combine_node[,2]=="TF")
pathway_combine_node=pathway_combine_node[c(tmp1,tmp2,tmp4,tmp3),]
pathway_combine_node=pathway_combine_node[!duplicated(pathway_combine_node[,1]),]

tmp1=pathway_combine_node[,1]
tmp2=rep("A",length(tmp1))
for (i in 1:length(tmp1)){
  for (j in 1:length(pathway_branch_final_node_attributes)){
    if (tmp1[i]%in%pathway_branch_final_node_attributes[[j]][,1])
      tmp2[i]=paste(tmp2[i],names(pathway_branch_final_node_attributes)[j],sep=";")
  }
  tmp2[i]=substr(tmp2[i],3,1000)
}
pathway_combine_node=cbind(pathway_combine_node,tmp2)
colnames(pathway_combine_node)[3]="pathway"

for (i in 1:length(pathway_branch_final_node_attributes)){
  pathway_branch_final_node_attributes[[i]]=merge(pathway_branch_final_node_attributes[[i]],pathway_combine_node[,c(1,3)],by.x=1,by.y=1,sort=F)
}

#save result

dir="./fisher_kegg"
dir_sub="./fisher_kegg/sub_pathway"
if(!file.exists(dir)){
  dir.create(dir)
  dir.create(dir_sub)
}else if(!file.exists(dir_sub)){
  dir.create(dir_sub)
}

pathway_overview=Pathway_fisher[filter_test]
pathway_overview=cbind(pathway_overview,paste(names(pathway_overview),"_nodes.txt",sep=""),paste(names(pathway_overview),".txt",sep=""))
pathway_overview[,1]=signif(as.numeric(pathway_overview[,1]),3)
write.table(pathway_overview,paste(dir,"/pathway_file_list.txt",sep=""),col.names=F,quote=F,sep="\t")


for (i in 1:length(pathway_branch_final_network)){
  write.table(pathway_branch_final_network[[i]],paste(dir,"/sub_pathway/",names(pathway_branch_final_network)[i],".txt",sep=""),row.names=F,quote=F,sep="\t")
  write.table(pathway_branch_final_node_attributes[[i]],paste(dir,"/sub_pathway/",names(pathway_branch_final_node_attributes)[i],"_nodes.txt",sep=""),row.names=F,quote=F,sep="\t")
}
write.table(pathway_combine,paste(dir,"/combined pathway.txt",sep=""),row.names=F,quote=F,sep="\t")
write.table(pathway_combine_node,paste(dir,"/combined pathway_nodes.txt",sep=""),row.names=F,quote=F,sep="\t")


temp1=nrow(pathway_combine_node)
temp2=nrow(pathway_combine)
write.table(paste('{"nodeCount":"',temp1,'","pathwayCount":"',temp2,'"}',sep=""),paste0(dir,"/report.dat"),row.names=F,col.names=F,quote=F,sep="\t")

