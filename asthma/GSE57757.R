setwd("D:\\R\\asthma\\GSE52742")

r1 <- getGEO("GSE57757",getGPL = F)
r2 <- r1[[1]]
r3 <- pData(r2)
r4 <- as.data.frame(exprs(r2))
head(r4)
gpl1 <- getGEO("GPL6246")
gpl2 <- Table(gpl1)
colnames(gpl2)
gpl3 <- gpl2[,c(1,10)]
gpl3$symbol <- str_split_fixed(gpl3$gene_assignment,pattern = "//",3)[,2]
ID <- rownames(r4)
r4$ID <- ID
eset1 <- merge(r4,gpl3,by="ID")
colnames(eset1)
class(eset1$symbol)
eset1$symbol <- as.character(eset1$symbol)
eset2 <- aggregate(eset1[,c(2:8)],by=list(eset1$symbol),FUN = mean)
rownames(eset2) <- eset2$Group.1
eset2 <- eset2[,-1]
head(eset2)
eset2 <- eset2[-1,]
boxplot(eset2)
eset3 <- as.data.frame(t(eset2))
colnames(eset3)
data <- as.data.frame(eset3$' Ccl6 ' )
colnames(data) <- "Ccl6"
data$group <- c(rep("ctl",4),rep("Allergic Asthma",3))
library(ggplot2)
library(ggpubr)
library(ggsignif)
p1 <- ggboxplot(data = data,x="group",y= "Ccl6",fill = "group",bxp.errorbar = T,
                bxp.errorbar.width = 0.5,xlab = "Disease group",
                ylab = "CCL6 expression (log2+1)",ggtheme = theme_light(),
                add = "jitter",palette = "cell")
p2 <- p1 + stat_compare_means()
p2
ggsave("CCL6 expression.pdf",width = 8,height = 6)
