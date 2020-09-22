#prepare

##read scMLnet result

scMLnet=read.csv("../covid_balf_sc/net/vis.node.csv")

Ligand=scMLnet$Node_name[scMLnet$node_type=="Lig"]
Receptor=scMLnet$Node_name[scMLnet$node_type=="Rec"]
TF=scMLnet$Node_name[scMLnet$node_type=="TF"]

##read KEGG pathway information (graphite package)

KEGG_all_edge=read.delim("database/KEGG_all_edge.txt")
KEGG_all_edge=as.matrix(KEGG_all_edge)
KEGG_pathway=unique(KEGG_all_edge[,5])
KEGG_edge_symbol=list()
for (i in 1:length(KEGG_pathway)){
  KEGG_edge_symbol[[i]]=KEGG_all_edge[KEGG_all_edge[,5]%in%KEGG_pathway[i],]
}
names(KEGG_edge_symbol)=KEGG_pathway

#main

##read activated pathway
pathway_select <- read.table("output/fisher_pathway.txt",header = T,sep = "\t")
pathway_select <- pathway_select$pathway[pathway_select$pval < 0.05]

##get node info in pathway

KEGG_node_symbol <- data.frame()
for (pw in pathway_select){
  
  x=KEGG_edge_symbol[[pw]]
  y=unique(c(x[,1],x[,2]))
  
  lig=intersect(Ligand,y)
  rec=intersect(Receptor,y)
  tf=intersect(TF,y)
  
  
  tmp <- data.frame(gene=c(lig,rec,tf),
                    type=c(rep("Lig",length(lig)),rep("Rec",length(rec)),rep("TF",length(tf))),
                    pathway=rep(pw,length(c(lig,rec,tf))))
  KEGG_node_symbol <- rbind(KEGG_node_symbol,tmp)

}
write.csv(KEGG_node_symbol,"pw_mlnet.csv")
