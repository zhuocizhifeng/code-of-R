setwd("D:\\R\\asthma\\GSE89809")
library(GEOquery)
library(stringr)
library(data.table)
library(ggplot2)
r1 <- getGEO("GSE89809",getGPL = F)
r2 <- r1[[1]]
r3 <- pData(r2)
r4 <- as.data.frame(exprs(r2))
gpl1 <- getGEO("GPL13158")
gpl2 <- Table(gpl1)
colnames(gpl2)
gpl3 <- gpl2[,c(1,11)]
ID <- rownames(r4)
eset1 <- r4[ID %in% gpl3$ID,] %>% cbind(gpl3)
colnames(eset1)
eset2 <- aggregate(eset1[1:145],by=list(eset1$`Gene Symbol`),FUN= mean)
rownames(eset2) <- eset2$Group.1
eset2 <- eset2[,-1]
boxplot(eset2)
library(limma)
eset3 <- as.data.frame(normalizeBetweenArrays(eset2))
boxplot(eset3)
table(r3$source_name_ch1)
#epithelial brushing
control_sample <- r3[r3$source_name_ch1 == "Epithelial brushings gene expression data from healthy control",]$geo_accession
Epi_sample1 <- r3[r3$source_name_ch1 == "Epithelial brushings gene expression data from moderate asthmatic",]$geo_accession
Epi_sample2 <- r3[r3$source_name_ch1 == "Epithelial brushings gene expression data from severe asthmatic",]$geo_accession
exp1 <- eset2[,control_sample]
exp2 <- eset2[,Epi_sample1]
exp3 <- eset2[,Epi_sample2]
data <- cbind(exp1,exp2)
data <- cbind(data,exp3)
grouplist <- c(rep("CTL",18),rep("AS",24))
data1 <- as.data.frame(t(data))
data1 <- as.data.frame(data1$CCL23)
colnames(data1) <- "CCL23"
data1$group <- grouplist
library(ggplot2)
library(ggpubr)
library(ggsignif)
p1 <- ggboxplot(data = data1,x="group",y= "CCL23",fill = "group",bxp.errorbar = T,
                bxp.errorbar.width = 0.5,xlab = "Disease group",
                ylab = "CCL23 expression (log2+1)",ggtheme = theme_light(),
                add = "jitter",palette = "lancet")
p2 <- p1 + stat_compare_means()
p2
ggsave("CCL23 expression.pdf",width = 8,height = 6)
design <- model.matrix(~0+factor(grouplist))
colnames(design) <- levels(factor(grouplist))
contrast.matrix <- makeContrasts("AS-CTL",levels = design)
fit <- lmFit(data,design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
DEG <- topTable(fit2,number = Inf,adjust.method = "fdr")
head(DEG,20)
DEG$gene <- rownames(DEG)
write.csv(DEG,"DEG of at term vs EOPE.csv")
diffgene <- DEG[with(DEG,(abs(logFC)>1 & adj.P.Val < 0.05)),]
head(diffgene)
write.csv(diffgene,file = "diffgene-at term vs EOPE.csv")