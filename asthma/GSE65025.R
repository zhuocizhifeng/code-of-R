setwd("D:\\R\\asthma\\GSE65025")
library(GEOquery)
library(stringr)
library(data.table)
r1 <- getGEO("GSE65204",getGPL = F)
r2 <- r1[[1]]
r3 <- pData(r2)
r4 <- as.data.frame(exprs(r2))
gpl1 <- getGEO("GPL14550")
gpl2 <- Table(gpl1)
colnames(gpl2)
gpl3 <- gpl2[,c(1,7)]
ID <- rownames(r4)
head(ID,20)
r4$ID <- ID
eset1 <- merge(r4,gpl3,by="ID")
colnames(eset1)
eset2 <- aggregate(eset1[2:70],by=list(eset1$GENE_SYMBOL),FUN= mean)
rownames(eset2) <- eset2$Group.1
eset2 <- eset2[,-1]
boxplot(eset2)
library(limma)
eset2 <- as.data.frame(normalizeBetweenArrays(eset2))
boxplot(eset2)
table(r3$`asthma:ch1`)
grouplist <- ifelse(r3$`asthma:ch1` == "TRUE","AS","CTL")
head(grouplist)
design <- model.matrix(~0+factor(grouplist))
colnames(design) <- levels(factor(grouplist))
contrast.matrix <- makeContrasts("AS-CTL",levels = design)
fit <- lmFit(eset2,design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
DEG <- topTable(fit2,number = Inf,adjust.method = "fdr")
head(DEG,20)
write.csv(DEG,"DEG of AS VS CTL.csv")
diffgene <- DEG[with(DEG,(abs(logFC)>1 & adj.P.Val < 0.05)),]
head(diffgene)
write.csv(diffgene,file = "diffgene of AS VS CTL.csv")
library(pheatmap)
choose_gene = head(rownames(diffgene),50)
choose_matrix = eset2[choose_gene,]
choose_matrix = t(scale(t(choose_matrix)))
annotation_col=data.frame(group=grouplist)
rownames(annotation_col)=colnames(eset2)
pdf("heatmap of top50 DEGs.pdf",height = 10,width = 10)
pheatmap(choose_matrix,color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         scale = "row",annotation_col=annotation_col,show_rownames = F,show_colnames = F,
         main = "TOP50 heatmap ",cluster_rows = T,cluster_cols = F)
dev.off()

library(ggplot2)
logFC_cutoff = 1
DEG$change = as.factor(ifelse(DEG$P.Value < 0.05 &
                                abs(DEG$logFC) > logFC_cutoff,
                              ifelse(DEG$logFC > logFC_cutoff,
                                     'UP','DOWN'),'NOT'))
this_title = paste0('Cutoff for logFC is ',round(logFC_cutoff,2),
                    '\nThe number of up gene is ',nrow(DEG[DEG$change == 'UP',]),
                    '\nThe number of down gene is ',nrow(DEG[DEG$change == 'DOWN',]))

g = ggplot(data = DEG,
           aes(x = logFC, y = -log10(P.Value),
               color = change)) +
  geom_point(alpha = 0.4, size = 1.75) +
  theme_set(theme_set(theme_bw(base_size = 10))) +
  xlab("log2 fold change") + ylab('-log10 p-value') +
  ggtitle(this_title) +theme(plot.title = element_text(size = 10,hjust = 0.5)) +
  scale_color_manual(values = c('blue','black','red'))
pdf("volcano of DEG.pdf",width = 10,height = 10)
print(g)
dev.off()


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
ee=enrichGO(gene = gene.df$ENTREZID,OrgDb = org.Hs.eg.db, pvalueCutoff =1, qvalueCutoff = 1, ont="all", readable =T)
GO=as.data.frame(ee)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
#保存富集结果
write.table(GO,file="GO_DEGs.txt",sep="\t",quote=F,row.names = F)
#定义显示Term数目
showNum=5
if(nrow(GO)<30){
  showNum=nrow(GO)
}
#柱状???
pdf(file="DEGs_GObarplot.pdf",width = 20,height =18)
bar=barplot(ee, drop = TRUE, showCategory =showNum,split="ONTOLOGY",color = colorSel)+ facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()
#气泡???
pdf(file="DEGs_GObubble.pdf",width = 12,height =10)
bub=dotplot(ee,showCategory = showNum, orderBy = "GeneRatio",split="ONTOLOGY", color = colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bub)
bub
dev.off()
#KEGG
library(DO.db)
install.packages("R.utils")
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
pdf(file="DEGs_KEGGbarplot.pdf",width = 18,height = 16)
barplot(kk, drop = TRUE, showCategory = 20, color = colorSel)
dev.off()