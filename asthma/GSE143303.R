setwd("D:\\R\\asthma\\GSE143303")
library(GEOquery)
library(stringr)
library(data.table)
r1 <- getGEO("GSE143303",getGPL = F)
r2 <- r1[[1]]
r3 <- pData(r2)
r4 <- as.data.frame(exprs(r2))
r2
r5 <- fread("GSE143303_non-normalized.txt")
colnames(r5)
eset1 <- as.data.frame(r5[,c(1,4,6,8,10,12,14,16,18,20,
                             22,24,26,28,30,32,34,36,38,40,
                             42,44,46,48,50,52,54,56,58,60,
                             62,64,66,68,70,72,74,76,78,80,
                             82,84,86,88,90,92,94,96,98,100,102,
                             104,106,108,110,112,114,116,118,120,122)])
head(eset1)
colnames(eset1)
colnames(eset1)[1] <- "ID"
eset1 <- merge(eset1,gpl3,by="ID")
colnames(eset1)
eset2 <- aggregate(eset1[2:61],by=list(eset1$Symbol),FUN = mean)
rownames(eset2) <- eset2$Group.1
eset2 <- eset2[,-1]
head(eset2)
eset3 <- log2(eset2+1)
colnames(eset3) <- colnames(r4)[1:60]
boxplot(eset3)
grouplist <- r3$`inflammatory phenotype:ch1`
library(limma)
eset3 <- as.data.frame(normalizeBetweenArrays(eset3))
library(WGCNA)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
data <- eset3
head(data)
colnames(data)
m.vars <- apply(data,1,var)
dat <- data[which(m.vars>quantile(m.vars,probs=seq(0,1,0.25))[4]),]
dim(dat)
head(dat)
dat <- dat[-c(1,2,3,4),]
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
traitData$Paucigranulocytic_asthma <- ifelse(grouplist == "Paucigranulocytic asthma",1,0)
traitData$Eosinophilic_asthma <- ifelse(grouplist == "Eosinophilic asthma",1,0)
traitData$Neutrophilic_asthma <- ifelse(grouplist == "Neutrophilic asthma",1,0)
traitData$Healthy_control <- ifelse(grouplist == "Healthy control",1,0)
datTraits <- traitData[,-1]
head(datTraits)
rownames(datTraits) <- rownames(datExpr0)
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
               cex.lab.x = 0.5,
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
      abline(v=0.5,h=0.2,col="red")
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
write.csv(geneInfo, file = "WGCNA_GS_MM.csv")
for (mod in 1:nrow(table(moduleColors)))
{  
  modules = names(table(moduleColors))[mod]
  probes = colnames(datExpr0)
  inModule = (moduleColors == modules)
  modGenes = probes[inModule]
  write.table(modGenes, file =paste0("WGCNA_",modules,".txt"),sep="\t",row.names=F,col.names=F,quote=F)
}
save.image(file = "WGCNA.rdata")
g1 <- geneInfo
head(g1)
table(g1$moduleColor)
# 重新计算模块 eigengenes
MEs = moduleEigengenes(datExpr0, moduleColors)$eigengenes
# 提取临床特征weight
weight = as.data.frame(datTraits$Neutrophilic_asthma);
names(weight) = "Neutrophilic_asthma"
# 在eigengenes模块中加入临床特征weight
MET = orderMEs(cbind(MEs, weight))
# 绘制eigengenes和临床特征weight之间的关系图
pdf("WGCNA_EigengeneNetworks.pdf",width = 12,height = 12)
plotEigengeneNetworks(MET, "", 
                      marDendro = c(0,4,1,2), 
                      marHeatmap = c(3,4,1,2), 
                      cex.lab = 0.5)
dev.off()
# 分别绘制                      
# 绘制树状图
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# 绘制热图
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)
#magenta red purple module
g2 <- g1[g1$moduleColor == "purple",]
g3 <- g2[abs(g2$GS.Neutrophilic_asthma)>0.2,]
g4 <- g3[abs(g3$MMpurple)>0.5,]
hubgene <- g4$probes
write.table(hubgene,file="WGCNA_hubgene of purple.txt",sep="\t",row.names=F,col.names=F,quote=F)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
colorSel <- "qvalue"
gene <- bitr(hubgene,fromType = "SYMBOL",toType = "ENTREZID",
             OrgDb = org.Hs.eg.db)
ee1 <- enrichGO(gene$ENTREZID,OrgDb = org.Hs.eg.db,
                ont = "all",pvalueCutoff = 0.05,readable = T)
go1 <- as.data.frame(ee1)
write.csv(go1,"GO of purple module hub gene.csv")
pdf(file="WGCNA_GObarplot of hubgene in purple.pdf",width = 10,height = 8)
bar=barplot(ee1,drop = TRUE, showCategory =5,split="ONTOLOGY",color = colorSel)+ facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()
R.utils::setOption("clusterProfiler.download.method","auto")
kk <- enrichKEGG(gene = gene$ENTREZID,
                 organism = "hsa",
                 pvalueCutoff = 0.05)
kk1 <- as.data.frame(kk)
write.csv(kk1,"KEGG of purple module hub gene.csv")
pdf(file="WGCNA_KEGGdotplot of hubgene in purple.pdf",width = 10,height = 8)
dotplot(kk,  orderBy = "GeneRatio", showCategory = 10, color = colorSel)
dev.off()

# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datExpr0, power = 6);
modules = "magenta"
probes = colnames(datExpr0)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule]
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);
chooseTopHubInEachModule(datExpr = datExpr0,moduleColors)


library(limma)
boxplot(eset3)
write.csv(eset3,"gene exp.csv")
head(grouplist)
group <- as.data.frame(grouplist)
colnames(group) <- "group"
group$sample <- colnames(eset3)
NEA_sample <- group[group$group == "Neutrophilic asthma",]$sample
control_sample <- group[group$group == "Healthy control",]$sample
exp1 <- eset3[,control_sample]
exp2 <- eset3[,NEA_sample]
eset4 <- cbind(exp1,exp2)
grouplist1 <- c(rep("CTL",13),rep("NEA",9))

design <- model.matrix(~0+factor(grouplist1))
colnames(design) <- levels(factor(grouplist1))
contrast.matrix <- makeContrasts("NEA-CTL",levels = design)
fit <- lmFit(eset4,design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
DEG <- topTable(fit2,number = Inf,adjust.method = "fdr")
head(DEG,20)
write.csv(DEG,"DEG of NEA vs CTL.csv")
diffgene <- DEG[with(DEG,(abs(logFC)>1 & adj.P.Val < 0.05)),]
head(diffgene)
write.csv(diffgene,file = "diffgene of UP VS N.csv")
library(pheatmap)
eset5 <- diffgene[order(diffgene$logFC,decreasing = T),]
choose_gene = head(rownames(eset4),50)
choose_matrix = eset4[choose_gene,]
choose_matrix = t(scale(t(choose_matrix)))
eset6 <- t(eset4)
eset6 <- as.data.frame(eset6)
annotation_col=data.frame(group=grouplist1)
rownames(annotation_col)=rownames(eset6)
p = pheatmap(choose_matrix,color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
             show_rownames = T,show_colnames = F,main = "Diffgene heatmap",
             annotation_col = annotation_col,cluster_rows = T,cluster_cols = F,
             annotation_legend = T,scale = "row")
library(ggplot2)
pdf("heatmap of diffgene.pdf",width = 8,height = 8)
print(p)
dev.off()

colnames(DEG)
logFC_cutoff = 1
DEG$change = as.factor(ifelse(DEG$adj.P.Val < 0.05 &
                                abs(DEG$logFC) > logFC_cutoff,
                              ifelse(DEG$logFC > logFC_cutoff,
                                     'UP','DOWN'),'NOT'))
this_title = paste0('Cutoff for logFC is ',round(logFC_cutoff,2),
                    '\nThe number of up gene is ',nrow(DEG[DEG$change == 'UP',]),
                    '\nThe number of down gene is ',nrow(DEG[DEG$change == 'DOWN',]))
library(ggplot2)
pdf("volcano of different expression of NEA VS CTL.pdf",width = 10,height = 8)
g = ggplot(data = DEG,
           aes(x = logFC, y = -log10(P.Value),
               color = change)) +
  geom_point(alpha = 0.4, size = 1.75) +
  theme_set(theme_set(theme_bw(base_size = 10))) +
  xlab("log2 fold change") + ylab('-log10 p-value') +
  ggtitle(this_title) +theme(plot.title = element_text(size = 10,hjust = 0.5)) +
  scale_color_manual(values = c("blue",'black','red')) + theme_light()
print(g)
dev.off()
colorSel <- "qvalue"
gene = rownames(diffgene)
library(clusterProfiler)
library(org.Hs.eg.db)
gene.df <- bitr(gene, fromType = 'SYMBOL',
                toType = c('ENSEMBL','ENTREZID'),
                OrgDb = org.Hs.eg.db)
head(gene.df)
ee=enrichGO(gene = gene.df$ENTREZID,OrgDb = org.Hs.eg.db,pvalueCutoff = 1, ont="all", readable =T)
GO=as.data.frame(ee)
GO=GO[(GO$pvalue<0.05 & GO$qvalue<0.05),]
write.csv(GO,file="GO of diffgene in UP VS N.csv")
showNum= 5
pdf(file="GObarplot of diffgene of up vs n.pdf",width = 12,height =10)
bar=barplot(ee,drop = TRUE, showCategory =showNum,split="ONTOLOGY",color = colorSel)+ facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()
pdf(file="GObubble of diffgene of up vs n.pdf",width = 18,height =16)
bub=dotplot(ee,showCategory = showNum, orderBy = "GeneRatio",split="ONTOLOGY", color = colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bub)
dev.off()

###KEGG enrichment analysis of DEGs
R.utils::setOption("clusterProfiler.download.method","auto")
kk <- enrichKEGG(gene = gene.df$ENTREZID,
                 organism = "hsa",
                 pvalueCutoff = 0.05)
KEGG=as.data.frame(kk)
KEGG
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(DEG$Gene[match(strsplit(x,"/")[[1]],as.character(DEG$))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<1 & KEGG$qvalue<1),]
write.table(KEGG,file="KEGG of diffgene in high and low CIP2A.txt",sep="\t",quote=F,row.names = F)
showNum=30
if(nrow(KEGG)<showNum){
  showNum=nrow(KEGG)
}
pdf(file="KEGGdotplot of diffgene up vs n.pdf",width = 10,height = 8)
dotplot(kk,  orderBy = "GeneRatio", showCategory = showNum, color = colorSel)
dev.off()

###GSEA GO and KEGG analysis of DEGs
DEG1 <- DEG[order(DEG$logFC,decreasing = T),]
genelist <- DEG1$logFC
names(genelist) <- rownames(DEG1)
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
pdf("GSEA GO dotplot.pdf",width=10,height = 10)
dotplot(gse.GO)
dev.off()
pdf("GSEA go GSEAplot.pdf",width = 10,height = 10)
p <- gseaplot2(gse.GO,1,pvalue_table = T)
print(p)
dev.off()
gse.kegg <- gseKEGG(genelist4,organism = "hsa",
                    pvalueCutoff = 0.05)
pdf("GSEA KEGG dotplot top10.pdf",width=10,height = 10)
dotplot(gse.kegg)
dev.off()
pdf("GSEA KEGG GSEAplot.pdf",width = 10,height = 10)
p <- gseaplot2(gse.kegg,1,pvalue_table = T)
print(p)
dev.off()


diffgene <- fread("diffgene of NEA VS CTL.csv")
purple <- fread("WGCNA_hubgene of purple.txt",header = F)
red <- fread("WGCNA_hubgene of red.txt",header = F)
magenta <- fread("WGCNA_hubgene of magenta.txt",header = F)
library(ggvenn)
library(ggplot2)
list1 <- list(DEG=diffgene$V1,purple_module=purple$V1)
list2 <- list(DEG=diffgene$V1,red_module=red$V1)
list3 <- list(DEG=diffgene$V1,magenta_module=magenta$V1)
ggvenn(list1,show_percentage = F,text_size = 5,fill_color = c("blue","purple"))
ggsave(filename = "DEG vs purple_module.pdf",width = 8,height = 6)
ggvenn(list2,show_percentage = F,text_size = 5,fill_color = c("blue","red"))
ggsave(filename = "DEG vs red_module.pdf",width = 8,height = 6)
ggvenn(list3,show_percentage = F,text_size = 5,fill_color = c("blue","magenta"))
ggsave(filename = "DEG vs magenta_module.pdf",width = 8,height = 6)

#NEA VS EOA
group <- as.data.frame(grouplist)
colnames(group) <- "group"
group$sample <- colnames(eset3)
NEA_sample <- group[group$group == "Neutrophilic asthma",]$sample
EOA_sample <- group[group$group == "Eosinophilic asthma",]$sample
exp1 <- eset3[,EOA_sample]
exp2 <- eset3[,NEA_sample]
eset4 <- cbind(exp1,exp2)
grouplist1 <- c(rep("EOA",22),rep("NEA",9))

design <- model.matrix(~0+factor(grouplist1))
colnames(design) <- levels(factor(grouplist1))
contrast.matrix <- makeContrasts("NEA-EOA",levels = design)
fit <- lmFit(eset4,design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
DEG <- topTable(fit2,number = Inf,adjust.method = "fdr")
head(DEG,20)
write.csv(DEG,"DEG of NEA vs EOA.csv")
diffgene <- DEG[with(DEG,(abs(logFC)>1 & adj.P.Val < 0.05)),]
head(diffgene)
write.csv(diffgene,file = "diffgene of NEA VS EOA.csv")

colnames(DEG)
logFC_cutoff = 1
DEG$change = as.factor(ifelse(DEG$adj.P.Val < 0.05 &
                                abs(DEG$logFC) > logFC_cutoff,
                              ifelse(DEG$logFC > logFC_cutoff,
                                     'UP','DOWN'),'NOT'))
this_title = paste0('Cutoff for logFC is ',round(logFC_cutoff,2),
                    '\nThe number of up gene is ',nrow(DEG[DEG$change == 'UP',]),
                    '\nThe number of down gene is ',nrow(DEG[DEG$change == 'DOWN',]))
library(ggplot2)
pdf("volcano of different expression of NEA VS EOA.pdf",width = 10,height = 8)
g = ggplot(data = DEG,
           aes(x = logFC, y = -log10(P.Value),
               color = change)) +
  geom_point(alpha = 0.4, size = 1.75) +
  theme_set(theme_set(theme_bw(base_size = 10))) +
  xlab("log2 fold change") + ylab('-log10 p-value') +
  ggtitle(this_title) +theme(plot.title = element_text(size = 10,hjust = 0.5)) +
  scale_color_manual(values = c("blue",'black','red')) + theme_light()
print(g)
dev.off()

#NEA VS PGA
group <- as.data.frame(grouplist)
colnames(group) <- "group"
group$sample <- colnames(eset3)
NEA_sample <- group[group$group == "Neutrophilic asthma",]$sample
PGA_sample <- group[group$group == "Paucigranulocytic asthma",]$sample
exp1 <- eset3[,PGA_sample]
exp2 <- eset3[,NEA_sample]
eset4 <- cbind(exp1,exp2)
grouplist1 <- c(rep("PGA",16),rep("NEA",9))

design <- model.matrix(~0+factor(grouplist1))
colnames(design) <- levels(factor(grouplist1))
contrast.matrix <- makeContrasts("NEA-PGA",levels = design)
fit <- lmFit(eset4,design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
DEG <- topTable(fit2,number = Inf,adjust.method = "fdr")
head(DEG,20)
write.csv(DEG,"DEG of NEA vs PGA.csv")
diffgene <- DEG[with(DEG,(abs(logFC)>1 & adj.P.Val < 0.05)),]
head(diffgene)
write.csv(diffgene,file = "diffgene of NEA VS PGA.csv")
library(ggplot2)
pdf("heatmap of diffgene.pdf",width = 8,height = 8)
print(p)
dev.off()

colnames(DEG)
logFC_cutoff = 1
DEG$change = as.factor(ifelse(DEG$adj.P.Val < 0.05 &
                                abs(DEG$logFC) > logFC_cutoff,
                              ifelse(DEG$logFC > logFC_cutoff,
                                     'UP','DOWN'),'NOT'))
this_title = paste0('Cutoff for logFC is ',round(logFC_cutoff,2),
                    '\nThe number of up gene is ',nrow(DEG[DEG$change == 'UP',]),
                    '\nThe number of down gene is ',nrow(DEG[DEG$change == 'DOWN',]))
library(ggplot2)
pdf("volcano of different expression of NEA VS PGA.pdf",width = 10,height = 8)
g = ggplot(data = DEG,
           aes(x = logFC, y = -log10(P.Value),
               color = change)) +
  geom_point(alpha = 0.4, size = 1.75) +
  theme_set(theme_set(theme_bw(base_size = 10))) +
  xlab("log2 fold change") + ylab('-log10 p-value') +
  ggtitle(this_title) +theme(plot.title = element_text(size = 10,hjust = 0.5)) +
  scale_color_manual(values = c("blue",'black','red')) + theme_light()
print(g)
dev.off()


ctl_gene <- fread("diffgene of NEA VS CTL.csv")$V1
EOA_gene <- fread("diffgene of NEA VS EOA.csv")$V1
PGA_gene <- fread("diffgene of NEA VS PGA.csv")$V1
list4 <- list(ctl=ctl_gene,eoa=EOA_gene,pga=PGA_gene,magenta=magenta$V1)
ggvenn(list4,show_elements = F,show_percentage = F)
ggsave(filename = "unique gene of magenta module among asthma group.pdf",width = 8,height = 6)

gene11 <- ctl_gene[ctl_gene %in% EOA_gene]
gene11 <- gene11[gene11 %in% PGA_gene]
gene22 <- gene11[gene11 %in% red$V1]
write.csv(file = "unique gene of red module among asthma group.csv",gene22)

gene33 <- gene11[gene11 %in% magenta$V1]
write.csv(file = "unique gene of magenta module among asthma group.csv",gene33)
colorSel <- "qvalue"
gene = gene33
library(clusterProfiler)
library(org.Hs.eg.db)
gene.df <- bitr(gene, fromType = 'SYMBOL',
                toType = c('ENSEMBL','ENTREZID'),
                OrgDb = org.Hs.eg.db)
head(gene.df)
ee=enrichGO(gene = gene.df$ENTREZID,OrgDb = org.Hs.eg.db,pvalueCutoff = 1, ont="all", readable =T)
GO=as.data.frame(ee)
GO=GO[(GO$pvalue<0.05 & GO$qvalue<0.05),]
write.csv(GO,file="GO of unique gene of magenta module.csv")
showNum= 5
pdf(file="GObarplot of unique gene of magenta module.pdf",width = 10,height =8)
bar=barplot(ee,drop = TRUE, showCategory =showNum,split="ONTOLOGY",color = colorSel)+ facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()

R.utils::setOption("clusterProfiler.download.method","auto")
kk <- enrichKEGG(gene = gene.df$ENTREZID,
                 organism = "hsa",
                 pvalueCutoff = 0.05)
KEGG=as.data.frame(kk)
KEGG
KEGG=as.data.frame(kk)
KEGG=KEGG[(KEGG$pvalue<1 & KEGG$qvalue<1),]
write.csv(KEGG,file="KEGG of unique gene of magenta module.csv")
showNum=10
pdf(file="KEGGdotplot of unique gene of magenta module.pdf",width = 10,height = 8)
dotplot(kk,  orderBy = "GeneRatio", showCategory = showNum, color = colorSel)
dev.off()


#logistic regression
gene33
head(eset3)
data1 <- as.data.frame(t(eset3))
data2 <- data1[,gene33]
data2$group <- grouplist
class(data2$group)
data2$group <- factor(data2$group,levels = c("Healthy control","Eosinophilic asthma",
                                             "Paucigranulocytic asthma","Neutrophilic asthma"),
                      labels = c("0","EA","PGA","NEA"))
summary(data2)
library(nnet)
form <- as.formula("group ~ .")
model <- multinom(form,data2,model = T,probabilities = T)
summary(model)
sink("logistic regression model.txt")
summary(model)
sink()
coef(model)
summary(model)
library(car)
sink("logistic regression model p_value.txt")
Anova(model)
sink()
#ADORA3 AQP9 CSF3R CCL20
library(questionr)
sink("logistic regression model OR.txt")
round(odds.ratio(model),2)
sink()

library(ggplot2)
library(ggpubr)
head(gene22)
data1 <- as.data.frame(t(eset3))
data2 <- data1[,gene22]
data2$group <- grouplist
class(data2$group)
data2$group <- factor(data2$group,levels = c("Healthy control","Eosinophilic asthma",
                                             "Paucigranulocytic asthma","Neutrophilic asthma"),
                      labels = c("CTL","EA","PGA","NEA"))
compare <- list(c("CTL","EA"),
                c("CTL","PGA"),
                c("CTL","NEA"),
                c("EA","PGA"),
                c("EA","NEA"),
                c("PGA","NEA"))

gene_colname <- colnames(data2)
for (i in 1:11) {
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

head(eset3)
NEA_sample
data_nea <- eset3[NEA_sample]
data_nea <- as.data.frame(t(data_nea))
data_plot <- data_nea$AIF1
data_plot <- as.data.frame(data_plot)
colnames(data_plot) <- "AIF1"
data_plot$IL2 <- data_nea$IL2
data_plot$IL10 <- data_nea$IL10
data_plot$IFN_gamma <- data_nea$IFNG
data_plot$TNF <- data_nea$TNF
data_11 <- as.data.frame(data_nea$AIF1)
colnames(data_11) <- "AIF1"
data_11$sample <- NEA_sample
data_11 <- data_11[order(data_11$AIF1,decreasing = F),]
data_11$sample <- factor(data_11$sample)
p1=ggplot(data_11,aes(x=sample,y=AIF1))
p2 <- p1 + geom_point(aes(color=AIF1,size=AIF1)) + labs(x="Patient ID",y="AIF1 expression level") +
  scale_color_gradient(low = "blue",high = "red")
library(pheatmap)
mycolors <- colorRampPalette(c("blue", "white", "red"), bias = 1.2)(100)
data_plot <- data_plot[order(data_plot$AIF1,decreasing = F),]
tmp=t(scale(data_plot))
tmp[tmp > 1] = 1
tmp[tmp < -1] = -1
p3=pheatmap(tmp,col= mycolors,show_colnames = F,cluster_cols = F,cluster_rows = F)
library(ggplotify)
p4 <- as.ggplot(as.grob(p3))
library(ggpubr)
ggarrange(p2,p4,heights = c(2,1),ncol = 1,nrow = 2)
plots = list(p2,p4)
lay1 = rbind(c(rep(1,7)),c(rep(2,7)))
library(gridExtra)
grid.arrange(grobs = plots, layout_matrix = lay1, heigths = c(2,2),weights=c(10,10))

head(data_plot)
data_plot$IL1B <- data_nea$NLRP3
data_plot$FCN1 <- data_nea$FCN1
data_plot$IL6<- data_nea$IL6
data_plot$TNFR1B <- data_nea$TNFRSF1B
p1 <- ggplot(data_plot,aes(x= TNFR1A,y=AIF1))
p1 + geom_point(size = 3,color = "red") + geom_smooth(method = "lm",color = "blue",fill = "lightgray") + 
  stat_cor(method = "pearson") + theme_light() + labs(x = "TNFR1A expression",y = "AIF1 expression")
ggsave("correlation plot of AIF1 and TNFR1A.pdf",width = 8,height = 6)

library(GSVA)
library(GSEABase)
library(limma)
library(stringr)
head(data_nea)
data_nea1 <- as.data.frame(t(data_nea))
geneset <- getGmt("h.all.v2022.1.Hs.symbols.gmt")
data1 <- as.matrix(data_nea1)
es <- gsva(data1,geneset,verbose=T,mx.diff=T)
rownames(es) <- str_replace(rownames(es),"HALLMARK_","")
es1 <- as.data.frame(t(es))
es1$AIF1 <- data_plot$AIF1
colnames(es1)
p1 <- ggplot(es1,aes(x= PANCREAS_BETA_CELLS,y=AIF1))
p1 + geom_point(size = 3,color = "red") + geom_smooth(method = "lm",color = "blue",fill = "lightgray") + 
  stat_cor(method = "pearson") + theme_light() + labs(x = "KRAS_SIGNALING_UP",y = "AIF1 expression")
colnames(es1)
ggsave("correlation plot of AIF1 and KRAS_SIGNALING_UP.pdf",width = 8,height = 6)


geneset <- getGmt("c2.all.v2022.1.Hs.symbols.gmt")
data1 <- as.matrix(data_nea1)
es <- gsva(data1,geneset,verbose=T,mx.diff=T)
es1 <- as.data.frame(t(es))
es1$AIF1 <- data_plot$AIF1
colnames(es1)
p1 <- ggplot(es1,aes(x= REACTOME_NF_KB_IS_ACTIVATED_AND_SIGNALS_SURVIVAL,y=AIF1))
p1 + geom_point(size = 3,color = "red") + geom_smooth(method = "lm",color = "blue",fill = "lightgray") + 
  stat_cor(method = "pearson") + theme_light() +
labs(x = "REACTOME_INNATE_IMMUNE_SYSTEM",y = "AIF1 expression")
ggsave("correlation plot of AIF1 and REACTOME_INNATE_IMMUNE_SYSTEM.pdf",width = 8,height = 6)


#normal
control_data <- eset3[control_sample]
control_data <- as.data.frame(t(control_data))
data_plot1 <- control_data$AIF1
data_plot1 <- as.data.frame(data_plot1)
colnames(data_plot1) <- "AIF1"
data_plot1$C1QC <- control_data$C1QC
data_plot1$NLRP3 <- control_data$NLRP3

p1 <- ggplot(data_plot1,aes(x= C1QC,y=AIF1))
p1 + geom_point(size = 3,color = "red") + geom_smooth(method = "lm",color = "blue",fill = "lightgray") + 
  stat_cor(method = "pearson") + theme_light() + labs(x = "NLRP3 expression",y = "AIF1 expression")
ggsave("correlation plot of AIF1 and TNFR1A.pdf",width = 8,height = 6)

library(GSVA)
library(GSEABase)
library(limma)
library(stringr)
head(data_nea1)
data_control <- as.data.frame(t(control_data))
head(data_control)
data1 <- cbind(data_nea1,data_control)
geneset <- getGmt("h.all.v2022.1.Hs.symbols.gmt")
data1 <- as.matrix(data1)
es <- gsva(data1,geneset,verbose=T,mx.diff=T)
exp4 <- as.data.frame(t(data1))
grouplist1 <- c(rep("NEA",9),rep("CTL",13))
design <- model.matrix(~0+factor(grouplist1))
colnames(design) <- levels(factor(grouplist1))
grouplist1
contrast.matrix <- makeContrasts("NEA-CTL",levels = design)
fit <- lmFit(es,design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
DEG2 <- topTable(fit2,number = Inf,adjust.method = "fdr")
head(DEG2,20)
write.csv(DEG2,"DEpathway(HALLMARK) of NEA VS CTL.csv")
diffgene2 <- DEG2[with(DEG2,(abs(t)>2 & adj.P.Val < 0.05)),]
head(diffgene2)
write.csv(diffgene2,file = "diff pathway(HALLMARK) of NEA VS CTL.csv")

dat_plot <- data.frame(id=rownames(DEG2),t=DEG2$t)
dat_plot$threshold = factor(ifelse(dat_plot$t  >-2, ifelse(dat_plot$t >= 2 ,'Up','NoSignifi'),'Down'),levels=c('Up','Down','NoSignifi'))
dat_plot <- dat_plot %>% arrange(t)
dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)
id <- str_replace(dat_plot$id,"HALLMARK_","")
dat_plot$id <- id
dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)
library(ggplot2)
library(RColorBrewer)
library(ggprism)
color <- brewer.pal(3,"Set1")
p <- ggplot(data = dat_plot,aes(x = id,y = t,fill = threshold)) +
  geom_col()+
  coord_flip() +
  scale_fill_manual(values = c('Up'= "#E41A1C",'NoSignifi'='#cccccc','Down'="#377EB8")) +
  geom_hline(yintercept = c(-2,2),color = 'white',size = 0.5,lty='dashed') +
  xlab('') +
  ylab('t value of GSVA score') + #注??????????转??
  guides(fill=F)+ # ????示图??
  theme_prism(border = T) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
p
low1 <- dat_plot %>% filter(t < -2) %>% nrow()
low0 <- dat_plot %>% filter( t < 0) %>% nrow()
high0 <- dat_plot %>% filter(t < 2) %>% nrow()
high1 <- nrow(dat_plot)
p <- p + geom_text(data = dat_plot[1:low1,],aes(x = id,y = 0.1,label = id),
                   hjust = 0,color = 'black') + # 小??-1??为??色??签
  geom_text(data = dat_plot[(low1 +1):low0,],aes(x = id,y = 0.1,label = id),
            hjust = 0,color = 'grey') + # ??色??签
  geom_text(data = dat_plot[(low0 + 1):high0,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'grey') + # ??色??签
  geom_text(data = dat_plot[(high0 +1):high1,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'black') # ????1??为??色??签
p
ggsave("gsva_pathway_hallmark gene_bar_ NEA VS CTL.pdf",p,width = 14,height  = 10)

head(es)
colnames(es)
head(diffgene2)
ID <- rownames(diffgene2)
es11 <- es[ID,]
es11 <- as.data.frame(t(es11))
normalize <- function(x){
  return((x-min(x))/(max(x)-min(x)))
}
es22 <- as.data.frame(lapply(es11, normalize))
es22$group <- grouplist1
data2 <- as.data.frame(t(data1))
es22$AIF1 <- data2$AIF1
path_name <- colnames(es22)
i=2
for (i in 1:16) {
  pathway_name <- paste0(path_name[i])
  filename1 <- paste0("Spearman correlation plot contain control and NEA with ",
                      pathway_name,".pdf")
  p1 <- ggplot(es22,aes_string(x= pathway_name,y= "AIF1"))
  p1 + geom_point(size = 3,aes(color = group)) + geom_smooth(method = "lm",color = "blue",fill = "lightgray") + 
    stat_cor(method = "spearman") + theme_light() +
    labs(x = pathway_name,y = "AIF1 expression") +
    scale_color_lancet()
  ggsave(filename = filename1,width = 8,height = 6)
}
p1 <- ggplot(es22,aes(x= HALLMARK_INTERFERON_ALPHA_RESPONSE,y= AIF1))
p1 + geom_point(size = 3,aes(color = group)) + geom_smooth(method = "lm",color = "blue",fill = "lightgray") + 
  stat_cor(method = "spearman") + theme_light() +
  labs(x = "REACTOME_INNATE_IMMUNE_SYSTEM",y = "AIF1 expression") +
  scale_color_lancet()
ggsave("correlation plot of AIF1 and REACTOME_INNATE_IMMUNE_SYSTEM.pdf",width = 8,height = 6)



es_na <- es11[c(1:9),]
colnames(es_na)
es_na <- es_na[,-17]
es_na <- as.data.frame(lapply(es_na, normalize))
es_na$AIF1 <- data_nea$AIF1
colnames(es_na)
i=1
pname <- colnames(es_na)[i]
p1 <- ggplot(es_na,aes(x= HALLMARK_ALLOGRAFT_REJECTION,y= AIF1))
p1 + geom_point(size = 3,color = "red") + geom_smooth(method = "lm",color = "blue",fill = "lightgray") + 
  stat_cor(method = "pearson") + theme_light() +
  labs(x = "REACTOME_INNATE_IMMUNE_SYSTEM",y = "AIF1 expression")
ggsave("correlation plot of AIF1 and REACTOME_INNATE_IMMUNE_SYSTEM.pdf",width = 8,height = 6)


