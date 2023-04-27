setwd("D:\\R\\asthma\\GSE45111")
library(GEOquery)
library(stringr)
library(data.table)
r1 <- getGEO("GSE45111",getGPL = F)
r2 <- r1[[1]]
r3 <- pData(r2)
r4 <- as.data.frame(exprs(r2))
r5 <- fread("GSE45111_non-normalized.txt",skip = 4)
sample<- r5[1,]
colnames(r5)
r5 <- r5[,c(1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,
            32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,
            62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94)]
r2
gpl1 <- getGEO("GPL6104")
gpl2 <- Table(gpl1)
colnames(gpl2)
gpl3 <- gpl2[,c(1,12)]
colnames(gpl3)[1] <- "ID_REF"
eset1 <- merge(r5,gpl3,by="ID_REF")
colnames(eset1)
eset2 <- aggregate(eset1[,2:48],by=list(eset1$Symbol),FUN = mean)
rownames(eset2) <- eset2$Group.1
eset2 <- eset2[,-1]
boxplot(eset2)
eset2 <- log2(eset2+1)
boxplot(eset2)
colnames(eset2) <- colnames(r4)
data <- as.data.frame(t(eset2))
magenta_gene <- c("SIGLEC10","CFB","SOD2","PDZK1IP1","UBD","TNFRSF6B",
                  "SAA1","ADORA3","AQP9","C1QA","FCGR1A",
                  "FGR","FPR1","CSF3R","DAPL1","CCL20")
data1 <- data[,c("SIGLEC10","CFB","SOD2","PDZK1IP1","UBD","TNFRSF6B",
                 "SAA1","ADORA3","AQP9","C1QA","FGR","FPR1","CSF3R","CCL20")]
data2 <- data[,c("LY96","IFI30","FKBP11","FCER1G","C1QB","C1QC","AIF1","TNFRSF1B",
                 "CD14","FCN1")]
library(ggplot2)
library(ggpubr)
data1$group <- r3$`asthma phenotype:ch1`
data1$group <- as.factor(data1$group)
compare <- list(c("Eosinophilic Asthma","Neutrophilic Asthma"),
                c("Eosinophilic Asthma","Paucigranulocytic Asthma"),
                c("Neutrophilic Asthma","Paucigranulocytic Asthma"))

p1 <- ggboxplot(data = data1,x="group",y= "CCL20",bxp.errorbar = T,
                bxp.errorbar.width = 0.2,xlab = "Asthma subgroup",
                ylab = "gene expression (log2+1)",ggtheme = theme_light(),
                add = "jitter",fill = "group",palette = "lancet")
p1 +stat_compare_means(comparisons = compare)
ggsave("CCL20 expression.pdf",width = 10,height = 8)

data2$group <- r3$`asthma phenotype:ch1`
data2$group <- as.factor(data2$group)
compare <- list(c("Eosinophilic Asthma","Neutrophilic Asthma"),
                c("Eosinophilic Asthma","Paucigranulocytic Asthma"),
                c("Neutrophilic Asthma","Paucigranulocytic Asthma"))

p1 <- ggboxplot(data = data2,x="group",y= "FCN1",bxp.errorbar = T,
                bxp.errorbar.width = 0.2,xlab = "Asthma subgroup",
                ylab = "gene expression (log2+1)",ggtheme = theme_light(),
                add = "jitter",fill = "group",palette = "lancet")
p1 +stat_compare_means(comparisons = compare)
ggsave("FCN1 expression.pdf",width = 10,height = 8)
