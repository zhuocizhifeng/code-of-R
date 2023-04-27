library(ggplot2)
library(GEOquery)
library(stringr)
library(data.table)
setwd("D:\\R\\asthma\\GSE22528")
r1 <- getGEO("GSE22528",getGPL = F)
r2 <- r1[[1]]
r3 <- pData(r2)
r4 <- as.data.frame(exprs(r2))
gpl1 <- getGEO("GPL96")
gpl2 <- Table(gpl1)
colnames(gpl2)
gpl3 <- gpl2[,c(1,11)]
ID <- rownames(r4)
eset1 <- r4[ID %in% gpl3$ID,] %>% cbind(gpl3)
colnames(eset1)
eset2 <- aggregate(eset1[1:10],by=list(eset1$`Gene Symbol`),FUN = mean)
rownames(eset2) <- eset2$Group.1
eset2 <- eset2[,-1]
head(eset2)
boxplot(eset2)
eset3 <- as.data.frame(t(eset2))
data <- as.data.frame(eset3$CCL23)
colnames(data) <- "CCL23"
data$group <- r3$`disease state:ch1`
library(ggplot2)
library(ggpubr)
library(ggsignif)
p1 <- ggboxplot(data = data,x="group",y= "CCL23",fill = "group",bxp.errorbar = T,
                bxp.errorbar.width = 0.5,xlab = "Disease group",
                ylab = "CCL23 expression (log2+1)",ggtheme = theme_light(),
                add = "jitter",palette = "cell")
p2 <- p1 + stat_compare_means()
p2
ggsave("CCL23 expression.pdf",width = 8,height = 6)
