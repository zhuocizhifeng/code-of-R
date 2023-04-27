library(GEOquery)
library(stringr)
library(data.table)
setwd("D:\\R\\asthma\\GSE74075")
r1 <- getGEO("GSE74075",getGPL = F)
r2 <- r1[[1]]
r2
r3 <- pData(r2)
r4 <- as.data.frame(exprs(r2))
head(r4)
r5 <- fread("GSE74075_non-normalized.txt")
colnames(r5)
r5 <- r5[-c(1:5),]

eset1 <- as.data.frame(r5[,c(1,seq(from = 2, to=32, by=2))])
head(eset1)
colnames(eset1)
colnames(eset1)[1] <- "ID"
gpl1 <- getGEO("GPL6883")
gpl2 <- Table(gpl1)
colnames(gpl2)
gpl3 <- gpl2[,c(1,12)]
eset1 <- merge(eset1,gpl3,by="ID")
colnames(eset1)
for (i in 2:17) {
  eset1[,i] <- as.numeric(eset1[,i])
  
}
eset2 <- aggregate(eset1[2:17],by=list(eset1$Symbol),FUN = mean)
rownames(eset2) <- eset2$Group.1
eset2 <- eset2[,-1]
head(eset2)
eset3 <- log2(eset2+1)
colnames(eset3) <- colnames(r4)[1:16]
boxplot(eset3)
grouplist <- r3$`disease state:ch1`
library(limma)
boxplot(eset3)
write.csv(eset3,"gene_exp.csv")
write.csv(r3,"patient data.csv")
design <- model.matrix(~0+factor(grouplist))
colnames(design) <- levels(factor(grouplist))
contrast.matrix <- makeContrasts("asthmatic-healthy",levels = design)
fit <- lmFit(eset3,design)
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
#GO for DEGs
library(clusterProfiler)
library(org.Hs.eg.db)
pvalueFilter=0.05
qvalueFilter=1
colorSel="qvalue"
#GO富集分析
gene <- rownames(diffgene)
#entrezIDs <- mget(gene, org.Hs.egSYMBOL2EG, ifnotfound=NA)
gene.df <- bitr(gene, fromType = 'SYMBOL',
                toType = c('ENSEMBL','ENTREZID'),
                OrgDb = org.Hs.eg.db)
head(gene.df)
ee=enrichGO(gene = gene.df$ENTREZID,OrgDb = org.Hs.eg.db,ont="all", readable =T)
GO=as.data.frame(ee)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
write.table(GO,file="GO_DEGs.txt",sep="\t",quote=F,row.names = F)

#柱状???
pdf(file="DEGs_GObarplot.pdf",width = 10,height =8)
bar=barplot(ee, drop = TRUE, showCategory = 5,split="ONTOLOGY",color = colorSel)+ facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()
#KEGG
library(DO.db)
R.utils::setOption("clusterProfiler.download.method","auto")
kk <- enrichKEGG(gene = gene.df$ENTREZID,
                 organism = "hsa",
                 pvalueCutoff = 0.05)
KEGG=as.data.frame(kk)
KEGG
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(DEG$Gene[match(strsplit(x,"/")[[1]],as.character(DEG$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
write.table(KEGG,file="KEGG.txt",sep="\t",quote=F,row.names = F)
pdf(file="DEGs_KEGGbarplot.pdf",width = 10,height = 8)
barplot(kk, drop = TRUE, showCategory = 10, color = colorSel)
dev.off()
