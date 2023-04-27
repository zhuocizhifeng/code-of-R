library(GEOquery)
library(stringr)
library(data.table)
setwd("D:\\R\\asthma\\GSE109455")
r1 <- getGEO("GSE109455",getGPL = F)
r2 <- r1[[1]]
r2
r3 <- pData(r2)
r4 <- as.data.frame(exprs(r2))
head(r4)
boxplot(r4[1:10])
r4 <- log2(r4+1)
r2
gpl1 <- getGEO("GPL10904")
gpl2 <- Table(gpl1)
colnames(gpl2)
rownames(r4)
colnames(r4)
eset1 <- aggregate(r4[1:196],by=list(rownames(r4)),FUN = mean)
rownames(eset1) <- eset1$Group.1
eset1 <- eset1[,-1]
boxplot(eset1[1:10])
celltype <- r3$source_name_ch1
CD14sample <- r3[r3$source_name_ch1 == "CD14+ leukocytes",]$geo_accession
CD4sample <- r3[!r3$source_name_ch1 == "CD14+ leukocytes",]$geo_accession
library(limma)
#CD14
CD14 <- r3[r3$source_name_ch1 == "CD14+ leukocytes",]
grouplist <- ifelse(CD4$`bmi:ch1` > 25,"obese","control")
design <- model.matrix(~0+factor(grouplist))
colnames(design) <- levels(factor(grouplist))
contrast.matrix <- makeContrasts("obese-control",levels = design)
eset2 <- eset1[CD14sample]
fit <- lmFit(eset2,design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
DEG <- topTable(fit2,number = Inf,adjust.method = "fdr")
head(DEG,20)
write.csv(DEG,"DEG of as vs ctl.csv")
diffgene <- DEG[with(DEG,(abs(logFC)>1 & adj.P.Val < 0.05)),]
head(diffgene)
write.csv(diffgene,file = "diffgene of as VS ctl.csv")
colnames(DEG)
DEG$change = as.factor(ifelse(DEG$adj.P.Val < 0.05 &
                                abs(DEG$logFC) > logFC_cutoff,
                              ifelse(DEG$logFC > logFC_cutoff,
                                     'UP','DOWN'),'NOT'))
logFC_cutoff = 1
this_title = paste0('Cutoff for logFC is ',round(logFC_cutoff,2),
                    '\nThe number of up gene is ',nrow(DEG[DEG$change == 'UP',]),
                    '\nThe number of down gene is ',nrow(DEG[DEG$change == 'DOWN',]))
library(ggplot2)
pdf("volcano of AS vs ctl.pdf",width = 14,height = 10)
g = ggplot(data = DEG,
           aes(x = logFC, y = -log10(P.Value),
               color = change)) +
  geom_point(alpha = 0.4, size = 1.75) +
  theme_set(theme_set(theme_bw(base_size = 10))) +
  xlab("log2 fold change") + ylab('-log10 p-value') +
  ggtitle(this_title) +theme(plot.title = element_text(size = 10,hjust = 0.5)) +
  scale_color_manual(values = c('blue','black','red'))
print(g)
dev.off()
head(diffgene)