library(stringr)
library(GEOquery)
library(ggplot2)
library(data.table)
setwd("D:\\R\\asthma\\GSE108417")
r1 <- getGEO("GSE108417",getGPL = F)
r2 <- r1[[1]]
r3 <- pData(r2)
r4 <- as.data.frame(exprs(r2))
boxplot(r4)
r2
gpl1 <- Table(getGEO("GPL6246"))
head(gpl1)
symbol <-str_split_fixed(gpl1$gene_assignment,pattern = "//",3)[,2]
symbol <- as.data.frame(symbol)
symbol$symbol <- str_replace(symbol$symbol," ","")
symbol$ID <- gpl1$ID
ID <- rownames(r4)
r4$ID <- ID
eset1 <- merge(r4,symbol,by = "ID")
head(eset1)
colnames(eset1)
eset2 <- aggregate(eset1[2:13],by=list(eset1$symbol),FUN = mean)
head(eset2)
rownames(eset2) <- eset2$Group.1
eset2 <- eset2[,-1]
eset2 <- eset2[-1,]
head(eset2)
boxplot(eset2)
library(limma)
eset3 <- normalizeBetweenArrays(eset2)
eset3 <- as.data.frame(eset3)
boxplot(eset3)
#ctl,ea,na
eset4 <- eset3[c(1,2,3,4,5,6,10,11,12)]
grouplist <- c("ctl","ctl","ctl","ea","ea","ea","na","na","na")
data <- as.data.frame(t(eset4))
colnames(data)
data1 <- as.data.frame(data$`Aif1 `)
colnames(data1) <- "Aif1"
data1$group <- grouplist
library(ggpubr)
compare <- list(c("ctl","ea"),c("ctl","na"),c("ea","na"))
p1 <- ggboxplot(data = data1,x="group",y= "Aif1",bxp.errorbar = T,
                bxp.errorbar.width = 0.2,xlab = "Aif1 group",
                ylab = "gene expression",ggtheme = theme_light(),
                add = "jitter",fill = "group",palette = "lancet")
p1 + stat_compare_means(comparisons = compare,method = "t.test")
ggsave(filename = "Aif1 expression.pdf",width = 8,height = 6)
