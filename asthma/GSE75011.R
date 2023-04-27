setwd("D:\\R\\asthma\\GSE75011")
library(WGCNA)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DESeq2)
library(data.table)
library(GEOquery)
rm(list = ls())
pdata1 <- getGEO("GSE75011",getGPL = F)
pdata1 <- pdata1[[1]]
pdata <- pData(pdata1)
##### input the expression table and data convert####
exp <- fread("GSE75011_Raw_counts.tsv") #download from UCSC Xena cohort: TCGA TARGET GTEx
exp <- as.data.frame.matrix(exp)
rownames(exp) <- exp[,1]
exp <- exp[,2:ncol(exp)]
pd <- fread("GSE75011_series_matrix.txt")
r1 <- getGEO("GSE75011",getGPL = F)
r2 <- r1[[1]]
r3 <- pData(r2)
table(r3$`disease group:ch1`)
grouplist <- c(rep("AS",40),rep("CTL",15))
name <- colnames(exp)
exp_c <- exp[,c(2,26:33,35:65,66:80)]
condition <- factor(grouplist,levels = c("AS","CTL"))
colData <- data.frame(row.names = colnames(exp_c),condition)
dds <- DESeqDataSetFromMatrix(countData = exp_c, colData = colData, design = ~ condition) #Generate the DESeq2 input file
dds <- dds[rowSums(counts(dds)) > 1,]#filter the low express genes
dds <- DESeq(dds) #differential analysis
deseq_normal_vst_expression <- assay(vst(dds,blind = T))
write.csv(deseq_normal_vst_expression,"gene exp of AS vs CTL after DESEQ2 normolized.csv")
res <- results(dds,contrast = c("condition","AS","CTL"),independentFiltering = F) #obtain the analysis result of DESeq2
resorder <- res[order(res$pvalue),]
DEG <- as.data.frame(resorder)
logfc_cutoff <- with(DEG,mean(abs(log2FoldChange))+2*sd(abs(log2FoldChange)))
logfc_cutoff
DEG$change <- as.factor(ifelse(DEG$padj< 0.05 & abs(DEG$log2FoldChange)> logfc_cutoff,
                               ifelse(DEG$log2FoldChange > logfc_cutoff,"UP","DOWN"),"NOT"))
write.csv(DEG,"DEG of AS vs CTL.csv")
diffgene <- DEG[with(DEG,(abs(DEG$log2FoldChange)>1 & DEG$padj < 0.05)),]
head(diffgene)
write.csv(diffgene,"diffgene of AS vs CTL.csv")
#heatmap
library(pheatmap)
eset1 <- deseq_normal_vst_expression
head(eset1)
colnames(eset1)
choose_gene = head(rownames(diffgene),50)
choose_matrix = eset1[choose_gene,]
choose_matrix = t(scale(t(choose_matrix)))
eset2 <- t(eset1[choose_gene,])
eset2 <- as.data.frame(eset2)
eset3 <- as.data.frame(t(eset1))

annotation_col=data.frame(group=grouplist,sample=colnames(eset1))
annotation_col1 <- annotation_col[order(annotation_col$group,decreasing = T),]
annotation_col2 <- as.data.frame(annotation_col1$group)
colnames(annotation_col2) <- "Group"
rownames(annotation_col2) <- annotation_col1$sample
eset2$group <- grouplist
eset2 <- eset2[order(eset2$group,decreasing = T),]
eset4 <- eset2[,-51]
eset5 <- t(eset4)
eset5 <- scale(eset5)
pdf("heatmap of TOP50.pdf",width = 8,height = 8)
p = pheatmap(eset5,color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
             show_rownames = F,show_colnames = F,main = "TOP50 heatmap",
             annotation_col = annotation_col2,cluster_rows = T,cluster_cols = F,
             annotation_legend = T)
dev.off()
#volcano
colnames(DEG)
logFC_cutoff = 1
this_title = paste0('Cutoff for logFC is ',round(logFC_cutoff,2),
                    '\nThe number of up gene is ',nrow(DEG[DEG$change == 'UP',]),
                    '\nThe number of down gene is ',nrow(DEG[DEG$change == 'DOWN',]))
library(ggplot2)
pdf("volcano of AS vs ctl.pdf",width = 14,height = 10)
g = ggplot(data = DEG,
           aes(x = log2FoldChange, y = -log10(pvalue),
               color = change)) +
  geom_point(alpha = 0.4, size = 1.75) +
  theme_set(theme_set(theme_bw(base_size = 10))) +
  xlab("log2 fold change") + ylab('-log10 p-value') +
  ggtitle(this_title) +theme(plot.title = element_text(size = 10,hjust = 0.5)) +
  scale_color_manual(values = c('blue','black','red'))
print(g)
dev.off()
#GO
library(ggplot2)
colorSel <- "qvalue"
gene = rownames(diffgene)
gene.df <- bitr(gene, fromType = 'SYMBOL',
                toType = c('ENSEMBL','ENTREZID'),
                OrgDb = org.Hs.eg.db)
head(gene.df)
ee=enrichGO(gene = gene.df$ENTREZID,OrgDb = org.Hs.eg.db,pvalueCutoff = 1, ont="all", readable =T)
GO=as.data.frame(ee)
GO=GO[(GO$pvalue<0.05 & GO$qvalue<0.05),]
write.table(GO,file="GO of diffgene in AS VS CTL.txt",sep="\t",quote=F,row.names = F)
showNum=5
pdf(file="GObarplot of diffgene in AS VS CTL.pdf",width = 10,height =10)
bar=barplot(ee,drop = TRUE, showCategory = 5,split="ONTOLOGY",color = colorSel)+ facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()

#KEGG
R.utils::setOption("clusterProfiler.download.method","auto")
kk <- enrichKEGG(gene = gene.df$ENTREZID,
                 organism = "hsa",
                 pvalueCutoff = 0.05)
KEGG=as.data.frame(kk)
KEGG
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(DEG$Gene[match(strsplit(x,"/")[[1]],as.character(DEG$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<1 & KEGG$qvalue<1),]
write.table(KEGG,file="KEGG of diffgene in AS VS CTL.txt",sep="\t",quote=F,row.names = F)
showNum=10
if(nrow(KEGG)<showNum){
  showNum=nrow(KEGG)
}
pdf(file="KEGGbarplot of diffgene in AS VS CTL.pdf",width = 10,height = 10)
barplot(kk,  orderBy = "GeneRatio", showCategory = showNum, color = colorSel)
dev.off()
#GSEA
eset1 <- DEG
eset2 <- eset1[order(eset1$log2FoldChange,decreasing = T),]
genelist <- eset2$log2FoldChange
names(genelist) <- rownames(eset2)
genelist1 <- as.data.frame(genelist)
genelist1$symbol <- rownames(genelist1)
gene1 <- bitr(genelist1$symbol,OrgDb = org.Hs.eg.db,fromType = "SYMBOL",toType = "ENTREZID")
colnames(gene1)[1] <- "symbol"
genelist2 <- merge(genelist1,gene1,by= "symbol")
genelist3 <- genelist2[order(genelist2$genelist,decreasing = T),]
genelist4 <- genelist3$genelist
names(genelist4) <- genelist3$ENTREZID
library(enrichplot)
gse.GO <- gseGO(
  genelist4,
  ont = "ALL",  
  OrgDb = org.Hs.eg.db, 
  keyType = "ENTREZID",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)
pdf("GSEA GO dotplot top10.pdf",width=10,height = 10)
dotplot(gse.GO)
dev.off()
pdf("GSEA go GSEAplot.pdf",width = 10,height = 10)
p <- gseaplot2(gse.GO,1:10)
print(p)
dev.off()
gse.kegg <- gseKEGG(genelist4,organism = "hsa",
                    pvalueCutoff = 1)
pdf("GSEA KEGG dotplot top10.pdf",width=10,height = 10)
dotplot(gse.kegg)
dev.off()
pdf("GSEA KEGG GSEAplot.pdf",width = 10,height = 10)
p <- gseaplot2(gse.kegg,1:5)
print(p)
dev.off()



#WGCNA
data <- deseq_normal_vst_expression
head(data)
colnames(data)
library(WGCNA)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
m.vars <- apply(data,1,var)
dat <- data[which(m.vars>quantile(m.vars,probs=seq(0,1,0.25))[4]),]
dim(dat)
head(dat)
datExpr0 <- t(dat)
gsg = goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = hclust(dist(datExpr0), method = "average")
traitData <- as.data.frame(grouplist)
table(traitData$grouplist)
traitData$AS <- ifelse(traitData$grouplist == "AS",1,0)
traitData$CTL <- ifelse(traitData$grouplist == "AS",0,1)
datTraits <- traitData
rownames(datTraits) <- colnames(data)
datTraits <- datTraits[,-1]
traitColors = numbers2colors(datTraits, signed = T)
sampleTree2 = hclust(dist(datExpr0), method = "average")
pdf("WGCNA_sample dendrogram and trait heatmap.pdf",width=12,height=12)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()
enableWGCNAThreads()   
powers = c(c(1:10), seq(from = 12, to=20, by=2))    
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
pdf(file="WGCNA_scale_independence.pdf",width=9,height=5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.85,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
ppower <- sft$powerEstimate
softPower = ppower 
cor <- WGCNA::cor
net = blockwiseModules(datExpr0, power = softPower,maxBlockSize = 100000,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "TOM",
                       verbose = 3)
cor <- stats::cor
table(net$colors)
sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
pdf(file="WGCNA_merged_dynamic.pdf", width = 9, height = 6)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "WGCNA1.RData")
nGenes = ncol(datExpr0);
nSamples = nrow(datExpr0);
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
sizeGrWindow(10,6)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
pdf(file="WGCNA_Module_trait.pdf",width=10,height=10)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = F,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.lab.y = 0.3,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
traitNames=names(datTraits)
geneTraitSignificance = as.data.frame(cor(datExpr0, datTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")
for (trait in traitNames){
  traitColumn=match(trait,traitNames)  
  for (module in modNames){
    column = match(module, modNames)
    moduleGenes = moduleColors==module
    if (nrow(geneModuleMembership[moduleGenes,]) > 1){
      outPdf=paste("WGCNA_", trait, "_", module,".pdf",sep="")
      pdf(file=outPdf,width=7,height=7)
      par(mfrow = c(1,1))
      verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                         abs(geneTraitSignificance[moduleGenes, traitColumn]),
                         xlab = paste("Module Membership in", module, "module"),
                         ylab = paste("Gene significance for ",trait),
                         main = paste("Module membership vs. gene significance\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
      abline(v=0.8,h=0.5,col="red")
      dev.off()
    }
  }
}
moduleColors = labels2colors(net$colors)
which.trait <- "weight_g"
y <- datTraits[, 1]
GS <- as.numeric(cor(y ,datExpr0, use="p"))
GeneSignificance <-  abs(GS)
ModuleSignificance <- tapply(
  GeneSignificance,
  moduleColors, mean, na.rm=T)
plotModuleSignificance(GeneSignificance, moduleColors)
probes = colnames(datExpr0)
geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)
for (Tra in 1:ncol(geneTraitSignificance))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneTraitSignificance[,Tra],
                         GSPvalue[, Tra])
  names(geneInfo0) = c(oldNames,names(geneTraitSignificance)[Tra],
                       names(GSPvalue)[Tra])
}

for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[,mod],
                         MMPvalue[, mod])
  names(geneInfo0) = c(oldNames,names(geneModuleMembership)[mod],
                       names(MMPvalue)[mod])
}
geneOrder =order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]
write.table(geneInfo, file = "wgcna_GS_MM.xls",sep="\t",row.names=F)
for (mod in 1:nrow(table(moduleColors)))
{  
  modules = names(table(moduleColors))[mod]
  probes = colnames(datExpr0)
  inModule = (moduleColors == modules)
  modGenes = probes[inModule]
  write.table(modGenes, file =paste0("WGCNA_",modules,".txt"),sep="\t",row.names=F,col.names=F,quote=F)
}
save.image(file = "WGCNA.rdata")
g1 <- read.csv("WGCNA_GS_MM.csv")
head(g1)
table(g1$moduleColor)
g2 <- g1[g1$moduleColor == "blue",]
g3 <- g2[abs(g2$GS.AS)>0.5,]
g4 <- g3[abs(g3$MMblue)>0.8,]
hubgene <- g4$probes
write.table(hubgene,file="WGCNA_hubgene of blue.txt",sep="\t",row.names=F,col.names=F,quote=F)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
colorSel <- "qvalue"
gene <- bitr(hubgene,fromType = "SYMBOL",toType = "ENTREZID",
             OrgDb = org.Hs.eg.db)
ee <- enrichGO(gene$ENTREZID,OrgDb = org.Hs.eg.db,
               ont = "all",pvalueCutoff = 0.05,readable = T)
go <- as.data.frame(ee)
write.csv(go,"WGCNA_blue module_GO analysis.csv")
pdf(file="WGCNA_GObarplot of hubgene in blue.pdf",width = 10,height =8)
bar=barplot(ee,drop = TRUE, showCategory = 5,split="ONTOLOGY",color = colorSel)+ facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()
R.utils::setOption("clusterProfiler.download.method","auto")
kk <- enrichKEGG(gene = gene$ENTREZID,
                 organism = "hsa",
                 pvalueCutoff = 0.05)
kk1 <- as.data.frame(kk)
write.csv(kk1,"WGCNA_blue module_KEGG analysis")
pdf(file="KEGGdotplot of hubgene in blue.pdf",width = 14,height = 12)
dotplot(kk,  orderBy = "GeneRatio", showCategory = 10, color = colorSel)
dev.off()
pdf(file="KEGGbarplot of hubgene in turquosie.pdf",width = 14,height = 12)
barplot(kk,  orderBy = "GeneRatio", showCategory = 10, color = colorSel)
dev.off()


