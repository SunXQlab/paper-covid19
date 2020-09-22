library(DESeq2)
library(limma)

options(stringsAsFactors = F)

##prepare
count <- read.table("./data/count.txt",header = T,sep = "\t")

count$median=apply(count[,-ncol(count)],1,median)
count=count[order(count$gene,count$median,decreasing = T),]
count=count[!duplicated(count$gene),]
count=count[,-9]

expr=count
expr=expr[!(rowSums(expr[,-ncol(expr)]==0) == 7),]

rownames(expr) <- expr$gene
expr=expr[,-ncol(expr)]
colnames(expr) <- c(paste0("N",1:3),paste(paste0("P",rep(1:2,each=2)),rep(1:2,2),sep = "_"))

##sam_info
group <- factor(c(rep("Normal",3),rep("Patient",4)))
group <- relevel(group, ref = "Normal")
pdata <- data.frame(group = group,
                    batch = factor(c(rep(1,3),rep(2:3,2))))
rownames(pdata) <- colnames(expr)

##main
expr=round(expr)
dds <- DESeqDataSetFromMatrix(countData = expr,
                              colData = pdata,
                              design = ~ group) 
dds <- DESeq(dds)

expr_norm <- counts(dds, normalized=TRUE)
boxplot(log1p(expr_norm))
write.table(expr_norm,file='data/norm.txt',col.names = T,quote = FALSE,sep='\t')

resultsNames(dds)
res <- results(dds)
summary(res)
res

resOrdered <- res[order(res$padj),]
resOrdered=as.data.frame(resOrdered)
resOrdered$gene <- rownames(resOrdered)
head(resOrdered)

write.table(resOrdered, './output/DE_genes.txt', row.names = F, quote = F, sep = '\t')
