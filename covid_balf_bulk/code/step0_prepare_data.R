options(stringsAsFactors = F)

samples <- list.files("./data/")
genes <- read.table(paste("./data",samples[1],sep = "/"),header = TRUE, sep = "\t")[,1]

count = matrix(ncol = length(samples),nrow = length(genes))
i = 0
for(sample in samples){
  print(sample)
  i = i+1
  tmp = read.table(paste("./data",sample,sep = "/"),header = TRUE, sep = "\t")
  count[,i] <- tmp[,6]
}
count <- as.data.frame(count)
colnames(count) = unlist(lapply(samples, function(sample){strsplit(sample,"[.]")[[1]][1]}))
count$gene <- genes
write.table(count,file='data/count.txt',col.names = T,quote = FALSE,sep='\t')

TPM = matrix(ncol = length(samples),nrow = length(genes))
i = 0
for(sample in samples){
  print(sample)
  i = i+1
  tmp = read.table(paste("./data",sample,sep = "/"),header = TRUE, sep = "\t")
  TPM[,i] <- tmp[,7]
}
TPM <- as.data.frame(TPM)
colnames(TPM) = unlist(lapply(samples, function(sample){strsplit(sample,"[.]")[[1]][1]}))
TPM$gene <- genes
write.table(TPM,file='data/TPM.txt',col.names = T,quote = FALSE,sep='\t')
