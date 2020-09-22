library(ggplot2)
library(ggsignif)
options(stringsAsFactors = F)

cor_rectf_scmlnet <- read.csv("./scMLnet/cor_RecTF.csv")
cor_rectf_scmlnet$software <- rep("scMLnet",nrow(cor_rectf_scmlnet))
cor_rectf_scmlnet$absR <- abs(cor_rectf_scmlnet$R)

cor_rectf_nichnet <- read.csv("./NichNet/cor_RecTF.csv")
cor_rectf_nichnet$software <- rep("NicheNet",nrow(cor_rectf_nichnet))
cor_rectf_nichnet$absR <- abs(cor_rectf_nichnet$R)

df_rectf <- rbind(cor_rectf_scmlnet,cor_rectf_nichnet)
p1 <- ggplot(df_rectf,aes(x=software,y=absR,fill=software)) +
  geom_violin(alpha=0.8,width=1) + 
  geom_boxplot(width=0.2,position=position_dodge(0.9)) +
  guides(fill=F) + theme_bw() +
  xlab(' ')+ylab('abs(PCC) of Receptor and TF') +
  geom_signif(comparisons = list(c("scMLnet", "NicheNet")),
              map_signif_level = TRUE, textsize=6) + 
  theme(
    axis.title.y = element_text(size = 16),
    axis.text = element_text(size = 14),
    panel.grid = element_blank()
  )


cor_TFTar_scmlnet <- read.csv("./scMLnet/cor_TFTar.csv")
cor_TFTar_scmlnet$software <- rep("scMLnet",nrow(cor_TFTar_scmlnet))
cor_TFTar_scmlnet$absR <- abs(cor_TFTar_scmlnet$R)

cor_TFTar_nichnet <- read.csv("./NichNet/cor_TFTar_.csv")
cor_TFTar_nichnet$software <- rep("NicheNet",nrow(cor_TFTar_nichnet))
cor_TFTar_nichnet$absR <- abs(cor_TFTar_nichnet$R)

df_tftar <- rbind(cor_TFTar_scmlnet,cor_TFTar_nichnet)
p2 <- ggplot(df_tftar,aes(x=software,y=absR,fill=software)) +
  geom_violin(alpha=0.8,width=1) + 
  geom_boxplot(width=0.2,position=position_dodge(0.9)) +
  guides(fill=F) + theme_bw() +
  xlab(' ')+ylab('abs(PCC) of TF and Target gene') +
  geom_signif(comparisons = list(c("scMLnet", "NicheNet")),
              map_signif_level = TRUE, textsize=6) + 
  theme(
    axis.title.y = element_text(size = 16),
    axis.text = element_text(size = 14),
    panel.grid = element_blank()
  )


library(ggpubr)
p3 <- ggarrange(p1,p2,ncol = 2,labels = c("A","B"),align = "hv",font.label = list(size=24))
dpi=300
tiff("./merge.tiff",width = dpi*13, height = dpi*6, units = "px",res = dpi)
print(p3)
dev.off()
