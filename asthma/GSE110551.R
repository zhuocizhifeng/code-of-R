library(GEOquery)
library(stringr)
library(data.table)
r1 <- getGEO("GSE110551",getGPL = F)
r2 <- r1[[1]]
r3 <- pData(r2)
r4 <- as.data.frame(exprs(r2))
head(r4)
library(hgu133plus2.db)
gpl1 <- toTable(hgu133plus2SYMBOL)
ID <- rownames(r4)
eset1 <- r4[ID %in% gpl1$probe_id,] %>% cbind(gpl1)
head(eset1)
colnames(eset1)
eset2 <- aggregate(eset1[1:156],by=list(eset1$symbol),FUN = mean)
rownames(eset2) <- eset2$Group.1
eset2 <- eset2[,-1]
head(eset2)
boxplot(eset2)
library(limma)
eset3 <- as.data.frame(normalizeBetweenArrays(eset2))
boxplot(eset3)
setwd("D:\\R\\asthma\\GSE110551")
write.csv(eset3,"gene_exp.csv")
write.csv(r3,"patient data.csv")

asthma_group <- r3[r3$`asthma status:ch1` == "Asthma",]
asthma_and_obese_sample <- asthma_group[asthma_group$`obesity status:ch1` == "Obese",]$geo_accession
asthma_and_non_obese_sample <- asthma_group[!asthma_group$`obesity status:ch1` == "Obese",]$geo_accession
healthy_group <- r3[!r3$`asthma status:ch1` == "Asthma",]
healthy_and_obese_sample <- healthy_group[healthy_group$`obesity status:ch1` == "Obese",]$geo_accession
healthy_and_non_obese_sample <- healthy_group[!healthy_group$`obesity status:ch1` == "Obese",]$geo_accession
#group: hno,ho,ano,ao
exp1 <- eset3[healthy_and_non_obese_sample]
exp2 <- eset3[healthy_and_obese_sample]
exp3 <- eset3[asthma_and_non_obese_sample]
exp4 <- eset3[asthma_and_obvs CTese_sample]
eset4 <- cbind(exp1,exp2)
eset4 <- cbind(eset4,exp3)
eset4 <- cbind(eset4,exp4)
grouplist <- c(rep("hno",39),rep("ho",39),rep("ano",39),rep("ao",39))
#ano and ao
design <- model.matrix(~0+factor(grouplist))
colnames(design) <- levels(factor(grouplist))
contrast.matrix <- makeContrasts("ano-ao",levels = design)
fit <- lmFit(eset4,design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
DEG <- topTable(fit2,number = Inf,adjust.method = "fdr")
head(DEG,20)
write.csv(DEG,"DEG of asthma patients with obese vs non_asthma.csv")
diffgene <- DEG[with(DEG,(abs(logFC)>1 & adj.P.Val < 0.05)),]
head(diffgene)
write.csv(diffgene,file = "diffgene of asthma patients with obese vs non_asthma.csv")
gene22 <- read.csv(file = "D:\\R\\asthma\\GSE74075\\diffgene of as VS ctl.csv")
gene22 <- gene22$X
head(gene22)
ID <- colnames(data2)
gene33 <- gene22[gene22 %in% ID]
data1 <- as.data.frame(t(eset4))
data2 <- data1[,gene33]
data2$group <- grouplist
class(data2$group)
data2$group <- factor(data2$group,levels = c("hno","ho",
                                             "ano","ao"),
                      labels = c("hno","ho","ano","ao"))
compare <- list(c("hno","ho"),
                c("hno","ano"),
                c("hno","ao"),
                c("ho","ano"),
                c("ho","ao"),
                c("ano","ao"))
gene_colname <- colnames(data2)
library(ggpubr)
library(ggplot2)
for (i in 1:30) {
  genename <- paste0(gene_colname[i])
  filename <- paste0(gene_colname[i],"_expression.pdf",collapse = "_")
  p1 <- ggboxplot(data = data2,x="group",y= genename ,bxp.errorbar = T,
                  bxp.errorbar.width = 0.2,xlab = "Asthma subgroup",
                  ylab = "gene expression (log2+1)",ggtheme = theme_light(),
                  add = "jitter",fill = "group",palette = "lancet")
  p1 +stat_compare_means(comparisons = compare,label = "p.signif")
  ggsave(filename = filename,width = 10,height = 8)
  
}
data2$CD24 <- data1$CD24
p1 <- ggboxplot(data = data2,x="group",y="CD24" ,bxp.errorbar = T,
                bxp.errorbar.width = 0.2,xlab = "Asthma subgroup",
                ylab = "gene expression (log2+1)",ggtheme = theme_light(),
                add = "jitter",fill = "group",palette = "lancet")
p1 +stat_compare_means(comparisons = compare,label = "p.signif")
ggsave(filename = filename,width = 10,height = 8)
data2$BEX5 <- data1$BEX5
p1 <- ggboxplot(data = data2,x="group",y="BEX5" ,bxp.errorbar = T,
                bxp.errorbar.width = 0.2,xlab = "Asthma subgroup",
                ylab = "gene expression (log2+1)",ggtheme = theme_light(),
                add = "jitter",fill = "group",palette = "lancet")
p1 +stat_compare_means(comparisons = compare,label = "p.signif")
ggsave(filename = filename,width = 10,height = 8)