options(stringsAsFactors = F)

gset1 <- read.csv("../DEG/geneset1.csv",header = F)
gset1 <- gset1[,-c(1:2)]
gset1 <- unlist(gset1)

gset2 <- read.csv("../DEG/geneset2.csv",header = F)
gset2 <- gset2[,-c(1:2)]
gset2 <- unlist(gset2)

gset <- intersect(gset1,gset2)
write.table(gset,"../DEG/geneset.txt",quote = F,sep = "\t",row.names = F,col.names = F)
