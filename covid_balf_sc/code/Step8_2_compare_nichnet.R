setwd("./compare")

##prepare
library(stringr)
options(stringsAsFactors = F)

expression = readRDS("../data/pat_data.Rdata")
barcode <- read.table("../data/pat.barcodetype.txt",header = T,sep = "\t")

exprs <- t(as.matrix(expression))

info <- data.frame(
  ifelse(barcode$Cluster == "Secretory",1,0),
  ifelse(barcode$Cluster != "Secretory",1,0),
  as.vector(barcode$Cluster),
  as.vector(barcode$Barcode),
  unlist(lapply(barcode$Barcode, function(x){ str_split(x,"_",simplify = T)[2] }))
)
colnames(info) <- c("classified  as receptor cell","classified  as sender cell",
                    "non-receptor cell type","cell","sample")
info$`non-receptor cell type`[info$`classified  as receptor cell` == 1] = 0

pat_expression <- list()
pat_expression$expression <- exprs
pat_expression$sample_info <- info

saveRDS(object = pat_expression, file = "./NichNet/data/pat_expression.rds")

##############################################################################################

library(nichenetr)
library(tidyverse)

options(stringsAsFactors = F)

##input

#the prior Ligand-target model
ligand_target_matrix = readRDS("./NichNet/dataset/ligand_target_matrix.rds")
ligand_target_matrix[1:5,1:5]

#Expression data
pat_expression = readRDS("./NichNet/data/pat_expression.rds")
expression = pat_expression$expression
sample_info = pat_expression$sample_info # contains meta-information about the cells

## Define expressed genes ##

sender_cell = unique(sample_info$`non-receptor cell type`)[-c(8,9)]
sender_ids = lapply(sender_cell, function(x){
  sample_info %>% filter(`non-receptor cell type` == x) %>% pull(cell)
})
names(sender_ids) = sender_cell
secretory_ids = sample_info %>% filter(`classified  as receptor cell` == 1) %>% pull(cell)

expressed_genes_secretory = expression[secretory_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 4] %>% names()
expressed_genes_sender = list()
for (sender in sender_cell) {
  #sender <- sender_cell[1]
  ids <- sender_ids[[sender]]
  expressed_genes_sender[[sender]] <- expression[ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 4] %>% names()
}
expressed_genes_sender = expressed_genes_sender %>% unlist() %>% unique()

# Check the number of expressed genes
length(expressed_genes_sender)
length(expressed_genes_secretory)

## Define the gene set and background gene ##

geneset =  read.table("./DEG/geneset.txt",header = F)
geneset <- unlist(geneset)
head(geneset)

background_expressed_genes = expressed_genes_secretory %>% .[. %in% rownames(ligand_target_matrix)]
head(background_expressed_genes)

## Define a set of potential ligands ##

lr_network = readRDS("./NichNet/dataset/lr_network.rds")

ligands = lr_network %>% pull(from) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_sender)

receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_secretory)

lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) 
head(lr_network_expressed)

potential_ligands = lr_network_expressed %>% pull(from) %>% unique()
head(potential_ligands)

## NicheNet’s ligand activity analysis ##

ligand_activities = predict_ligand_activities(geneset = geneset, 
                                              background_expressed_genes = background_expressed_genes, 
                                              ligand_target_matrix = ligand_target_matrix, 
                                              potential_ligands = potential_ligands)

ligand_activities %>% arrange(-pearson) 

#show histogram of ligand activity scores
p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) +
  geom_histogram(color="black", fill="darkorange")  +
  # geom_density(alpha=.1, fill="orange") +
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(20, pearson) %>% pull(pearson))), color="red", linetype="dashed", size=1) +
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()
p_hist_lig_activity

best_upstream_ligands = ligand_activities %>% top_n(34, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
head(best_upstream_ligands)

## Infer and visualize target genes of top-ranked ligands ##

# Get the active ligand-target links
active_ligand_target_links_df = best_upstream_ligands %>% 
  lapply(get_weighted_ligand_target_links,
         geneset = geneset, 
         ligand_target_matrix = ligand_target_matrix, 
         n = 250) %>% 
  bind_rows()

##inferring the signaling paths

## input ##

weighted_networks = readRDS("./NichNet/dataset/weighted_networks.rds")
ligand_tf_matrix = readRDS("./NichNet/dataset/ligand_tf_matrix.rds")

lr_network = readRDS("./NichNet/dataset/lr_network.rds")
sig_network = readRDS("./NichNet/dataset/signaling_network.rds")
gr_network = readRDS("./NichNet/dataset/gr_network.rds")


ligand_target_links <-  active_ligand_target_links_df #read.table("./NichNet/output/new/active_ligand_target_links.txt",header = T,sep = "\t")

targets <- unique(ligand_target_links$target)

ligand_target_net <- list()
df_LigTar <- c()
for (target in targets) {
  
  #target <-targets[1]
  Lig <- unique(ligand_target_links$ligand[ligand_target_links$target == target])
  ligand_target_net[[target]] <- Lig
  df_LigTar <- c(df_LigTar,paste(Lig,target,sep = "_"))
  
}

## main ##

act_sig_net <- list()
for (i in 1:length(targets)) {
  
  print(i)
  target <-targets[i]
  ligand <- ligand_target_net[[target]]
  
  active_signaling_network = get_ligand_signaling_path(ligand_tf_matrix = ligand_tf_matrix,
                                                       ligands_all = ligand,
                                                       targets_all = target,
                                                       weighted_networks = weighted_networks)
  act_sig_net[[target]] <- active_signaling_network
  
}

## get the Rec-TF, TF-Tar pairs in paths

library(igraph)

Rec = lr_network %>% select(to) %>% unique()
TF = gr_network %>% select(from) %>% unique()

rec_tf_tar <- list()
for (i in 1:length(targets)){
  
  target <- targets[i]
  act_sig_net_tmp <- act_sig_net[[target]]
  
  g <- graph_from_data_frame(act_sig_net_tmp$sig, directed=TRUE)
  d  <- shortest.paths(g,mode="out",weights = NA)
  
  rec <- act_sig_net_tmp$sig %>% select(to) %>% intersect(Rec) %>% unlist() %>% unique()
  tf <- act_sig_net_tmp$gr %>% select(from) %>% intersect(TF) %>% unlist() %>% unique()
  
  if(length(rec) <= 0 | length(tf) <= 0) next
  
  is.link <- c()
  for (i in rec) {
    for (j in tf) {
      is.link <- c(is.link,!is.infinite(d[i,j]))
    }
  }
  names(is.link) <- paste(rep(rec,each=length(tf)),tf,sep = "_")
  
  sig_Rec_tf <- names(is.link)[is.link]
  if(length(sig_Rec_tf)>0){
    rec_tf_tar[[target]] <- sig_Rec_tf
  }
}

#####################################################################################################

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
network <- rec_tf_tar

#get cor of RecTF 
df_RecTF <- unlist(lapply(network, function(x){
  #x=network[[1]]
  rec=x$rec
  tf=x$tf
  paste(rec,tf,sep = "_")
}))

cor_RecTF <- getCorR(GCMat,df_RecTF,cells)
write.csv(cor_RecTF,file = "./compare/NichNet/cor_RecTF.csv")

#get cor of TFTar
df_TFTar <- c()
for (name in names(network)) {
  #name = "AGT"
  tf=network[[name]]$tf
  df_TFTar <- c(df_TFTar,paste(tf,name,sep = "_"))
}

cor_TFTar <- getCorR(GCMat,df_TFTar,cells)
write.csv(cor_TFTar,file = "./compare/NichNet/cor_TFTar.csv")

