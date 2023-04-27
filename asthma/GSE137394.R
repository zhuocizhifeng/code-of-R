setwd("D:\\R\\asthma\\GSE137394")
library(GEOquery)
library(stringr)
library(data.table)
r1 <- getGEO("GSE137394",getGPL = F)
r2 <- r1[[1]]
r3 <- pData(r2)
r4 <- as.data.frame(exprs(r2))
table(r3$`disease:ch1`)
library(hgu133plus2.db)
gpl1 <- toTable(hgu133plus2SYMBOL)
ID <- rownames(r4)
r4$probe_id <- ID
eset1 <- merge(r4,gpl1,by="probe_id")
colnames(eset1)
eset2 <- aggregate(eset1[,2:310],by=list(eset1$symbol),FUN= mean)
rownames(eset2) <- eset2$Group.1
eset2 <- eset2[,-1]
group <- as.data.frame(r3$`blood_eosinophil_count:ch1`)
group <- na.omit(group)
group$sample <- r3$geo_accession
group <- na.omit(group)
colnames(group)[1] <- "Count"
class(group$Count)
group$Count <- as.numeric(group$Count)
group <- na.omit(group)
group$status <- ifelse(group$Count < 500 , "N","UP")
table(group$status)
normal_sample <- group[group$status == "N",]$sample
up_sample <- group[group$status == "UP",]$sample
exp1 <- eset2[,normal_sample]
exp2 <- eset2[,up_sample]
eset3 <- cbind(exp1,exp2)
grouplist <- c(rep("N",238),rep("UP",58))
library(limma)
boxplot(eset3)
eset3 <- as.data.frame(normalizeBetweenArrays(eset3))
boxplot(eset3)
write.csv(eset3,"gene exp of UP vs N.csv")
head(grouplist)
design <- model.matrix(~0+factor(grouplist))
colnames(design) <- levels(factor(grouplist))
contrast.matrix <- makeContrasts("UP-N",levels = design)
fit <- lmFit(eset3,design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
DEG <- topTable(fit2,number = Inf,adjust.method = "fdr")
head(DEG,20)
write.csv(DEG,"DEG of UP vs N.csv")
diffgene <- DEG[with(DEG,(abs(logFC)>1 & adj.P.Val < 0.05)),]
head(diffgene)
write.csv(diffgene,file = "diffgene of UP VS N.csv")
library(pheatmap)
eset1 <- as.data.frame(deseq_normal_vst_expression)
head(eset1)
colnames(eset1)
eset4 <- diffgene[order(diffgene$logFC,decreasing = T),]
choose_gene = head(rownames(eset4),11)
choose_matrix = eset3[choose_gene,]
choose_matrix = t(scale(t(choose_matrix)))
eset5 <- t(eset3)
eset5 <- as.data.frame(eset5)
annotation_col=data.frame(group=grouplist)
rownames(annotation_col)=rownames(eset5)
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
pdf("volcano of different expression of UP vs N.pdf",width = 10,height = 8)
g = ggplot(data = DEG,
           aes(x = logFC, y = -log10(P.Value),
               color = change)) +
  geom_point(alpha = 0.4, size = 1.75) +
  theme_set(theme_set(theme_bw(base_size = 10))) +
  xlab("log2 fold change") + ylab('-log10 p-value') +
  ggtitle(this_title) +theme(plot.title = element_text(size = 10,hjust = 0.5)) +
  scale_color_manual(values = c('black','red'))
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
traitData$Normal <- ifelse(grouplist == "N",1,0)
traitData$UP <- ifelse(grouplist == "UP",1,0)
datTraits <- traitData[,-1]
head(datTraits)
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
g2 <- g1[g1$moduleColor == "salmon",]
g3 <- g2[abs(g2$GS.UP)>0.2,]
g4 <- g3[abs(g3$MMsalmon)>0.5,]
hubgene <- g4$probes
write.table(hubgene,file="WGCNA_hubgene of salmon.txt",sep="\t",row.names=F,col.names=F,quote=F)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
colorSel <- "qvalue"
gene <- bitr(hubgene,fromType = "SYMBOL",toType = "ENTREZID",
             OrgDb = org.Hs.eg.db)
ee1 <- enrichGO(gene$ENTREZID,OrgDb = org.Hs.eg.db,
               ont = "all",pvalueCutoff = 1,readable = T)
go1 <- as.data.frame(ee1)
pdf(file="WGCNA_GObarplot of hubgene in salmon.pdf",width = 10,height = 8)
bar=barplot(ee1,drop = TRUE, showCategory =5,split="ONTOLOGY",color = colorSel)+ facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()
R.utils::setOption("clusterProfiler.download.method","auto")
kk <- enrichKEGG(gene = gene$ENTREZID,
                 organism = "hsa",
                 pvalueCutoff = 0.05)
kk1 <- as.data.frame(kk)
pdf(file="KEGGdotplot of hubgene in turquosie.pdf",width = 14,height = 12)
dotplot(kk,  orderBy = "GeneRatio", showCategory = 10, color = colorSel)
dev.off()
pdf(file="KEGGbarplot of hubgene in turquosie.pdf",width = 14,height = 12)
barplot(kk,  orderBy = "GeneRatio", showCategory = 10, color = colorSel)
dev.off()

gene1 <- read.csv("diffgene of UP VS N.csv")
gene1 <- gene1$X
gene2 <- hubgene
list1 <- list(DEG=gene1,salmon=gene2)
library(ggvenn)
pdf("venn plot of DEGs and salmon.pdf",width = 8,height = 8)
ggvenn(list1,show_elements = F,show_percentage = F,fill_color = c("red","salmon"))
dev.off()

#immune cell of PBMC
library(stringr)
library(GSVA)
library(dplyr)
library(tibble)
geneSet <- read.csv("immune cell.csv") #downloaded from TISIDB
class(geneSet)
list <- as.list(geneSet)
a <- geneSet
a <- a[1:nrow(a),]
set <- colnames(a)
l <- list()
for (i in set) {
  x <-  as.character(a[,i])
  x <- x[nchar(x)!=0]
  x <-  as.character(x)
  l[[i]] <-x
}
library(limma)
r1 <- data
r2 <- r1
r2 <- as.matrix(r2)
class(list)
ssgsea <- gsva(r2,l,method="ssgsea",kcdf="Gaussian",abs.ranking=T)
ssgsea.1 <- ssgsea
for (i in colnames(ssgsea)) {
  
  ssgsea.1[,i] <- (ssgsea[,i] -min(ssgsea[,i]))/(max(ssgsea[,i] )-min(ssgsea[,i] ))
  
}
apply(ssgsea.1[,1:6], 2, range)
library(pheatmap)
pheatmap(ssgsea.1)
names <- colnames(ssgsea)
ssgsea.2 <- ssgsea
for (i in colnames(ssgsea)) {
  
  ssgsea.2[,i] <- (ssgsea[,i]/sum(ssgsea[,i])*100)
  
}
ssgsea.2 <- t(ssgsea.2)
class(ssgsea.2)
ssgsea.2 <- as.data.frame(ssgsea.2)
ssgsea.2$sample <- rownames(ssgsea.2)
library(ggsci)
library(tidyr)
library(ggpubr)
ssgsea.3 <- t(ssgsea.1)
ssgsea.3 <- as.data.frame(ssgsea.3)
grouplist <- group1$grouplist
colnames(ssgsea.3) <- group1$group
annotion_col=c(rep("low",257),rep("high",257))
pdf("immune heatmap.pdf",height = 9,width = 7)
pheatmap(ssgsea.1,cluster_cols = F,show_colnames = F,cluster_rows = F,
         annotation_col = grouplist)
dev.off()
library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
dat <- ssgsea.2 %>% as.data.frame() %>%
  rownames_to_column("Sample") %>%
  gather(key = Cell_type,value = Proportion,-Sample)
dat <- dat[-c(7105:7400),]
dat$Proportion <- as.numeric(dat$Proportion)
pdf(file="immune cell of N vs UP.pdf",width = 12,height = 7)
ggplot(dat,aes(Sample,Proportion,fill = Cell_type)) +
  geom_bar(stat = "identity") +
  labs(fill = "Cell Type",x ="Sample",y = "Immune cell proportion") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(28))
dev.off()

ggplot(dat,aes(Cell_type,Proportion,fill = Cell_type)) +
  geom_boxplot(outlier.shape = 21,color = "black") +
  theme_bw() +
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_fill_manual(values = mypalette(28))
dat$Group = grouplist
library(ggpubr)
pdf(file="difference plot of immune cell in N vs up.pdf",width = 9,height = 7)
ggplot(dat,aes(Cell_type,Proportion,fill = Group)) +
  geom_boxplot(outlier.shape = 21,color = "black") +
  theme_bw() +
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") +
  theme(axis.text.x = element_text(angle=90,vjust = 0.55))+
  scale_fill_manual(values = mypalette(28)[c(6,1)])+ stat_compare_means(aes(group = Group,label = ..p.signif..),method = "wilcox.test")
dev.off()
library(ggcorrplot)
pgsea <- cor_pmat(t(ssgsea))
corgsea <- round(cor(t(ssgsea)),3)
pdf("cor_heatmap_total.pdf",width = 10,height = 10)
ggcorrplot(corgsea,hc.order = F,  
           ggtheme = ggplot2::theme_void(base_size = 15), 
           colors = c("CornflowerBlue","White","tomato"), 
           lab = F,lab_size = 3,    
           tl.cex = 10,             
           p.mat = pgsea,        
           sig.level = 0.05,
           insig = "blank")       
dev.off()
sig_gene <- gene1
immu_data <- t(ssgsea)
datat <- t(data)
datat <- as.data.frame(datat)
library(psych)
x <- datat[,sig_gene]
y <- immu_data
library(psych)
d <- corr.test(x,y,use="complete",method = 'spearman')
r <- d$r
p <- d$p
library(ggcorrplot)
pdf("cor_heatmap_total gene vs immune.pdf",width = 14,height = 10)
ggcorrplot(t(d$r), show.legend = T,
           p.mat = t(d$p), digits = 2,  sig.level = 0.05,insig = "blank",lab = T,
           lab_size = 3)
dev.off()
library(ggcorrplot)
color <- brewer.pal(8,"Set1")
data1 <- data[sig_gene,]
pgsea <- cor_pmat(t(data1))
corgsea <- round(cor(t(data1)),3)
pdf("cor_heatmap_total_top10gene.pdf",width = 10,height = 10)
ggcorrplot(corgsea,hc.order = F,  
           ggtheme = ggplot2::theme_void(base_size = 15), 
           colors = c("CornflowerBlue","White","tomato"), 
           lab = T,lab_size = 3,    
           tl.cex = 10,            
           p.mat = pgsea,         
           sig.level = 0.05,
           insig = "blank")            
dev.off()


library(ggcorrplot)
ssgsea1 <- ssgsea[,normal_sample]
pgsea <- cor_pmat(t(ssgsea1))
corgsea <- round(cor(t(ssgsea1)),3)
pdf("cor_heatmap_normal.pdf",width = 10,height = 10)
ggcorrplot(corgsea,hc.order = F,  
           ggtheme = ggplot2::theme_void(base_size = 15), 
           colors = c("CornflowerBlue","White","tomato"), 
           lab = F,lab_size = 3,    
           tl.cex = 10,             
           p.mat = pgsea,        
           sig.level = 0.05,
           insig = "blank")       
dev.off()
sig_gene <- gene1
immu_data <- t(ssgsea1)
datat <- t(exp1)
datat <- as.data.frame(datat)
library(psych)
x <- datat[,sig_gene]
y <- immu_data
library(psych)
d <- corr.test(x,y,use="complete",method = 'spearman')
r <- d$r
p <- d$p
library(ggcorrplot)
pdf("cor_heatmap_eso_normal gene vs immune.pdf",width = 14,height = 10)
ggcorrplot(t(d$r), show.legend = T,
           p.mat = t(d$p), digits = 2,  sig.level = 0.05,insig = "blank",lab = T,
           lab_size = 3)
dev.off()
library(ggcorrplot)
color <- brewer.pal(8,"Set1")
data1 <- exp1[sig_gene,]
pgsea <- cor_pmat(t(data1))
corgsea <- round(cor(t(data1)),3)
pdf("cor_heatmap_eos_normal_gene.pdf",width = 10,height = 10)
ggcorrplot(corgsea,hc.order = F,  
           ggtheme = ggplot2::theme_void(base_size = 15), 
           colors = c("CornflowerBlue","White","tomato"), 
           lab = T,lab_size = 3,    
           tl.cex = 10,            
           p.mat = pgsea,         
           sig.level = 0.05,
           insig = "blank")            
dev.off()

#boxplot

head(eset3)
head(grouplist)
b_data <- as.data.frame(t(eset3))
b_data <- as.data.frame(b_data$CCL23)
colnames(b_data) <- "CCL23"
b_data$group <- grouplist
colnames(eset3)
c_data <- as.data.frame(t(eset2))
c_data <- as.data.frame(c_data$CCL23)
colnames(c_data) <- "CCL23"
c_data$EOS.count <- r3$`blood_eosinophil_count:ch1`
c_data <- na.omit(c_data)
c_data$EOS.count <- as.numeric(c_data$EOS.count)
c_data <- na.omit(c_data)
c_data$EOS.count <- c_data$EOS.count
data_t <- as.data.frame(t(eset3))
b_data$IL5 <- data_t$IL5
b_data$IL13 <- data_t$IL13
b_data$IL4 <- data_t$IL4
b_data$VCAM1 <- data_t$VCAM1
b_data$ICAM1 <- data_t$ICAM1
b_data$IL3 <- data_t$IL3
b_data$TNF <- data_t$TNF
b_data$IL17 <- data_t$IL17A
library(ggplot2)
library(ggpubr)
library(ggsignif)
p1 <- ggboxplot(data = b_data,x="group",y= "CCL23",fill = "group",bxp.errorbar = T,
          bxp.errorbar.width = 0.2,xlab = "eosinophila count",
          ylab = "CCL23 expression (log2+1)",ggtheme = theme_light(),
          add = "jitter",palette = "lancet")
p2 <- p1 + stat_compare_means()
p2
ggsave(filename = "CCL23 expression boxplot.pdf",plot = p2,width = 8,height = 6)
p3 <- ggscatter(data = c_data,x="EOS.count",y= "CCL23",add = "reg.line",
                color = "black",size=2,
                add.params = list(color= "blue",fill= "lightgray"),
                conf.int = T,cor.coef = T,
                cor.coeff.args = list(method = "spearman",
                                      label.x = 1000, label.sep = "\n",label.y = 15),
          xlab = "Eosinophil count",ylab = "CCL23 expression ",palette = "lancet")
p3 + theme_light()
ggsave(filename = "CCL23 and Eosinophil correlation.pdf",width = 10,height = 8)



p3 <- ggscatter(data = b_data,x="IL17",y= "CCL23",add = "reg.line",
                color = "black",size=2,
                add.params = list(color= "blue",fill= "lightgray"),
                conf.int = T,cor.coef = T,
                cor.coeff.args = list(method = "spearman",
                                      label.x = 3.5, label.sep = "\n",label.y = 12),
                xlab = "IL17A expression",ylab = "CCL23 expression ",palette = "cell")
p3 + theme_light()
ggsave(filename = "CCL23 and IL17A expression correlation.pdf",width = 10,height = 8)


head(b_data)
d_data <- b_data
d_data$group <- ifelse(d_data$CCL23 < median(d_data$CCL23),"low","high")
d_data$IL1B <- data_t$IL1B
d_data$IL9 <- data_t$IL9
d_data$IL6 <- data_t$IL6
d_data$IL18 <- data_t$IL18
p1 <- ggboxplot(data = d_data,x="group",y= "IL5",fill = "group",bxp.errorbar = T,
                bxp.errorbar.width = 0.2,xlab = "CCL23 expression",
                ylab = "IL5 expression (log2+1)",ggtheme = theme_light(),
                add = "jitter",palette = "lancet")
p2 <- p1 + stat_compare_means()
p2
ggsave(filename = "IL5 expression boxplot in CCL23 grouping.pdf",plot = p2,width = 8,height = 6)


new_grouplist <- d_data$group
design <- model.matrix(~0+factor(new_grouplist))
colnames(design) <- levels(factor(new_grouplist))
contrast.matrix <- makeContrasts("high-low",levels = design)
fit <- lmFit(eset3,design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
DEG <- topTable(fit2,number = Inf,adjust.method = "fdr")
head(DEG,20)
write.csv(DEG,"DEG of high ccl23 vs low ccl23.csv")
diffgene <- DEG[with(DEG,(abs(logFC)>1 & adj.P.Val < 0.05)),]
head(diffgene)
write.csv(diffgene,file = "diffgene of high ccl23 vs low ccl23.csv")
library(pheatmap)
eset1 <- as.data.frame(deseq_normal_vst_expression)
head(eset1)
colnames(eset1)
eset4 <- diffgene[order(diffgene$logFC,decreasing = T),]
choose_gene = head(rownames(eset4),11)
choose_matrix = eset3[choose_gene,]
choose_matrix = t(scale(t(choose_matrix)))
eset5 <- t(eset3)
eset5 <- as.data.frame(eset5)
annotation_col=data.frame(group=grouplist)
rownames(annotation_col)=rownames(eset5)
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
pdf("volcano of high ccl23 vs low ccl23.pdf",width = 10,height = 8)
g = ggplot(data = DEG,
           aes(x = logFC, y = -log10(P.Value),
               color = change)) +
  geom_point(alpha = 0.4, size = 1.75) +
  theme_set(theme_set(theme_bw(base_size = 10))) +
  xlab("log2 fold change") + ylab('-log10 p-value') +
  ggtitle(this_title) +theme(plot.title = element_text(size = 10,hjust = 0.5)) +
  scale_color_manual(values = c('black','red'))
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
