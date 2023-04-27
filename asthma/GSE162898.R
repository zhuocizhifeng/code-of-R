setwd("D:\\R\\asthma\\GSE162898")
library(hdf5r)
library(Seurat)
library(tidyverse)
data1 <- Read10X_h5("GSE162898_WT_filtered_feature_bc_matrix.h5")
data2 <- data1[[1]]
data3 <- data1[[2]]
head(data3)
data4 <- as.data.frame(data3)
data4 <- as.data.frame(t(data4))
r2 <- CreateSeuratObject(data2, project = "HDM_sample")
head(r2)
r2$Hto_lungmacrophage <- data4$HTO_lungMacrophage 
r2$blood_monocyte <- data4$HTO_bloodMonocyte
r2[["percent.mt"]] <- PercentageFeatureSet(r2,pattern = "^mt-")
head(r2)
dim(r2)
VlnPlot(r2,features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3)
r2 <- subset(r2,subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
dim(r2)

r3 <- NormalizeData(r2,normalization.method = "LogNormalize",scale.factor = 10000)
head(r3)
r3 <- FindVariableFeatures(r3,selection.method = "vst",nfeatures = 2000)
top10 <- head(VariableFeatures(r3),10)
plot3 <- VariableFeaturePlot(r3)
plot4 <- LabelPoints(plot = plot3,points = top10,repel = T)
plot3 + plot4
r3 <- ScaleData(r3,features = rownames(r3))
r3 <- RunPCA(r3,features = VariableFeatures(object = r3))
print(r3[["pca"]],dims = 1:5,nfeatures = 5)
VizDimLoadings(r3,dims = 1:2,reduction = "pca")
DimPlot(r3,reduction = "pca")
DimHeatmap(r3,dims = 1,cells = 500,balanced = T)
DimHeatmap(r3,dims = 1:15,cells = 500,balanced = T)
r3 <- JackStraw(r3,num.replicate = 100)
r3 <- ScoreJackStraw(r3,dims = 1:20)
JackStrawPlot(r3,dims = 1:20)
ElbowPlot(r3)
r3 <- FindNeighbors(r3,dims = 1:20)
r3 <- FindClusters(r3,resolution = 0.5)
r3 <- RunUMAP(r3,dims = 1:20)
head(r3)
table(r3@meta.data[["orig.ident"]])
DimPlot(r3,reduction = "umap",pt.size = 0.8,split.by = "orig.ident")
markers <- FindAllMarkers(r3,only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)
markers %>% group_by(cluster) %>% top_n(n = 2,wt = avg_log2FC)
top10m <- markers %>% group_by(cluster) %>% top_n(n = 10,wt = avg_log2FC)
DoHeatmap(r3,features = top10m$gene) + NoLegend()
top10m[top10m$cluster == 1,]$gene
VlnPlot(r3,features = c("Gatm"),slot = "counts",log = T)
VlnPlot(r3,features = top10m[top10m$cluster == 5,]$gene,slot = "counts",log = T)
FeaturePlot(r3,features = c("Aif1"),
            cols = c("gray","red")) 
DotPlot(r3,features = top10m[top10m$cluster == 1,]$gene)
DotPlot(r3,features = "Aif1")
#0:AM spp1,gpnmb,Lpl,Car4
#1:AM Car4,plet1,Lpl,Ear2
#2:IM C1qb,H2-DMb1
#3:Monocyte
#4:Monocyte
#5:Monocyte
#6:Macrophage
#7:Macrophage
#8:IM
#9:IM
#10:Monocyte
#11:AM siglecf
#12:IM
#13:IM
#14:Monocyte
#15:?
#16:IM  
#17:Monocyte Pf4
#18:Monocyte Alox15 Prg4
#19:Macrophage
#20:Monocyte ly6c1
#cell marker database
#new.cluster.ids1 <- c("C1_Progenitor osteoclasts","C2_Mature osteoclasts","C3_Dysfunctional osteoclasts")
#names(new.cluster.ids1) <- levels(r5)
#r5 <- RenameIdents(r5,new.cluster.ids1)
#DimPlot(r5,reduction = "umap",label = F)
library(SingleR)
library(celldex)
meta <- r3@meta.data
immune_mus_gene <- celldex::ImmGenData()
exp_singleR<- GetAssayData(r3, slot="data") 
m2.hpca <- SingleR(test = exp_singleR,ref = immune_mus_gene,labels = immune_mus_gene$label.main)
head(m2.hpca)
r3@meta.data$labels <- m2.hpca$labels
head(r3)
DimPlot(r3,reduction = "umap",pt.size = 0.5)
r3@meta.data$cluseter <- Idents(r3)
head(r3)
head(r3@meta.data)
meta <- r3@meta.data
DimPlot(r3, group.by = c("labels"),reduction = "umap",split.by = "orig.ident",pt.size = 0.8)

r4 <- subset(r3,subset = labels == c("Macrophages","Monocytes"))
head(r4)
meta1 <- r4@meta.data
DimPlot(r4,group.by = c("labels"),reduction = "umap",split.by = "orig.ident",pt.size = 0.8)
ggsave("Dimplot of mono and macro.pdf",width = 8,height = 6)
FeaturePlot(r4,features = c("Aif1"),
            cols = c("gray","red"),pt.size = 0.8) 
ggsave("Aif1 expression in mono and macro.pdf",width = 8,height = 6)
exp_singleR<- GetAssayData(r4, slot="data") 
m2.hpca <- SingleR(test = exp_singleR,ref = immune_mus_gene,labels = immune_mus_gene$label.fine)
head(m2.hpca)
Idents(r4)
DimPlot(r4,reduction = "umap")
rM <- subset(r4,subset = labels == "Macrophages")
head(rM)
DimPlot(rM)
rM <- NormalizeData(rM,normalization.method = "LogNormalize",scale.factor = 10000)
head(rM)
rM <- FindVariableFeatures(rM,selection.method = "vst",nfeatures = 2000)
top10 <- head(VariableFeatures(rM),10)
plot3 <- VariableFeaturePlot(r3)
plot4 <- LabelPoints(plot = plot3,points = top10,repel = T)
plot3 + plot4
rM <- ScaleData(rM,features = rownames(rM))
rM <- RunPCA(rM,features = VariableFeatures(object = rM))
print(rM[["pca"]],dims = 1:5,nfeatures = 5)
VizDimLoadings(rM,dims = 1:2,reduction = "pca")
DimPlot(rM,reduction = "pca")
DimHeatmap(rM,dims = 1,cells = 500,balanced = T)
DimHeatmap(rM,dims = 1:15,cells = 500,balanced = T)
rM <- JackStraw(rM,num.replicate = 100)
rM <- ScoreJackStraw(rM,dims = 1:20)
JackStrawPlot(rM,dims = 1:20)
ElbowPlot(rM)
rM <- FindNeighbors(rM,dims = 1:20)
rM <- FindClusters(rM,resolution = 0.5)
rM <- RunUMAP(rM,dims = 1:20)
head(rM)
table(rM@meta.data[["orig.ident"]])
DimPlot(rM,reduction = "umap",pt.size = 0.8)
M_markers <- FindAllMarkers(rM,only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)
M_markers %>% group_by(cluster) %>% top_n(n = 2,wt = avg_log2FC)
M_top10m <- M_markers %>% group_by(cluster) %>% top_n(n = 10,wt = avg_log2FC)
DoHeatmap(r3,features = M_top10m$gene) + NoLegend()

VlnPlot(rM,features = c("Itgax"),slot = "counts",log = T)
VlnPlot(rM,features = M_top10m[M_top10m$cluster == 0,]$gene,slot = "counts",log = T)
FeaturePlot(rM,features = c("Nos2"),
            cols = c("gray","red"),pt.size = 0.8) 
DotPlot(rM,features = M_top10m[M_top10m$cluster == 1,]$gene)
#Cd11c+ IM,Cd11- IM,AM
#0,1,2,3 AM (based on Itgam(Cd11b expression)),4,5,6,7,8,9,10
#Cd11c high expression 0,1,5,7,8,9
new.cluster.ids <- c("Cd11c+ IM","Cd11c+ IM","Cd11c- IM",
                     "AM","Cd11c- IM","Cd11c+ IM",
                     "Cd11c- IM","Cd11c+ IM","Cd11c+ IM","Cd11c+ IM",
                     "Cd11c- IM")
names(new.cluster.ids) <- levels(rM)
rM <- RenameIdents(rM,new.cluster.ids)
DimPlot(rM,reduction = "umap",pt.size = 0.8,label = F)
ggsave("Dimplot of Macrophage.pdf",width = 8,height = 6)
FeaturePlot(rM,features = c("Aif1"),
            cols = c("grey","red"),label = T,label.size = 5,label.color = "black") 
ggsave("Aif1 expression in macrophagy.pdf",width = 8,height = 6)
VlnPlot(rM,features = c("Aif1"),slot = "counts",log = T)
ggsave("Violin plor of aif1 expression in macrophage.pdf",width = 8,height = 6)
M_markers.1 <- FindAllMarkers(rM,only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)
M_top10m.1 <- M_markers.1 %>% group_by(cluster) %>% top_n(n = 10,wt = avg_log2FC)
DotPlot(rM,features = M_top10m.1[M_top10m.1$cluster == "Cd11c- IM",]$gene)

#GSVA
library(GSVA)
library(GSEABase)
gmt <- getGmt("mh.all.v2022.1.Mm.symbols.gmt")
data1 <- as.matrix(rM[["RNA"]]@data)
gsva1 <- gsva(data1,gmt,kcdf = "Gaussian")
head(gsva1)
head(rownames(gsva1))
pheatmap::pheatmap(gsva1,show_rownames = 1,show_colnames = 0)
gsva2 <- data.frame(t(gsva1),stringsAsFactors = F)
rM <- AddMetaData(rM,gsva2)
head(rM)
meat2 <- rM@meta.data
FeaturePlot(rM,features = c("Ccr2","Aif1"),reduction = "umap",
            cols = c("grey","blue","red"),label = T,pt.size = 0.8)
rownames(data1)
data2 <- data.frame(t(data1),stringsAsFactors = F)
colnames(data2)
data3 <- data2$Aif1
data3 <- as.data.frame(data3)
colnames(data3) <- "Aif1"
data3$Ccr2 <- data2$Ccr2
data3$group <- Idents(rM)
rownames(data3) <- rownames(data2)
data4 <- data3[data3$group == "Cd11c- IM",]
data4$Aif1group <- ifelse(data4$Aif1 < median(data4$Aif1),"low","high")
cd11cNsample <- rownames(data4)
cd11cNdata <- data1[,cd11cNsample]
gsva_cd11cN <- gsva(cd11cNdata,gmt,kcdf = "Gaussian")
library(limma)
grouplist1 <- data4$Aif1group
design <- model.matrix(~0+factor(grouplist1))
colnames(design) <- levels(factor(grouplist1))
grouplist1
contrast.matrix <- makeContrasts("high-low",levels = design)
fit <- lmFit(gsva_cd11cN,design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
DEG2 <- topTable(fit2,number = Inf,adjust.method = "fdr")
head(DEG2,20)
write.csv(DEG2,"DEpathway(HALLMARK) of aif1 high vs low in cd11c- IM.csv")
diffgene2 <- DEG2[with(DEG2,(abs(t)>2 & adj.P.Val < 0.05)),]
head(diffgene2)
write.csv(diffgene2,file = "diff pathway(HALLMARK) of aif1 high VS low in cd11c- IM.csv")

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
  ylab('t value of GSVA score') + #ע??????????ת??
  guides(fill=F)+ # ????ʾͼ??
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
                   hjust = 0,color = 'black') + # С??-1??Ϊ??ɫ??ǩ
  geom_text(data = dat_plot[(low1 +1):low0,],aes(x = id,y = 0.1,label = id),
            hjust = 0,color = 'grey') + # ??ɫ??ǩ
  geom_text(data = dat_plot[(low0 + 1):high0,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'grey') + # ??ɫ??ǩ
  geom_text(data = dat_plot[(high0 +1):high1,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'black') # ????1??Ϊ??ɫ??ǩ
p
ggsave("gsva_pathway_hallmark gene_bar_aif1 high VS low in cd11c- IM.pdf",p,width = 14,height  = 10)

color1 <- brewer.pal(4,"Set1")
col <- pal_lancet(palette = "lanonc")(2)

ID <- rownames(diffgene2)
es11 <- gsva_cd11cN[ID,]
es11 <- as.data.frame(t(es11))
normalize <- function(x){
  return((x-min(x))/(max(x)-min(x)))
}
es22 <- as.data.frame(lapply(es11, normalize))
es22$group <- grouplist1
es22$group <- as.factor(es22$group)
es22$Aif1 <- data4$Aif1
path_name <- colnames(es22)
i=1
for (i in 1:41) {
  pathway_name <- paste0(path_name[i])
  filename1 <- paste0("Spearman correlation plot of CD11c- IM contain AIF1 with ",
                      pathway_name,".pdf")
  p1 <- ggplot(es22,aes_string(x= pathway_name,y= "Aif1"))
  p1 + geom_point(size = 3,aes(color = group)) + geom_smooth(method = "lm",color = "blue",fill = "lightgray") + 
    stat_cor(method = "spearman") + theme_light() +
    labs(x = pathway_name,y = "AIF1 expression") +
    scale_color_manual(values = c("#ED0000FF","#00468BFF"))
  ggsave(filename = filename1,width = 8,height = 6)
}
p1 <- ggplot(es22,aes(x= HALLMARK_INTERFERON_ALPHA_RESPONSE,y= AIF1))
p1 + geom_point(size = 3,aes(color = group)) + geom_smooth(method = "lm",color = "blue",fill = "lightgray") + 
  stat_cor(method = "spearman") + theme_light() +
  labs(x = "REACTOME_INNATE_IMMUNE_SYSTEM",y = "AIF1 expression") +
  scale_color_lancet()
ggsave("correlation plot of AIF1 and REACTOME_INNATE_IMMUNE_SYSTEM.pdf",width = 8,height = 6)



rM0 <- subset(r4,subset = labels == "Monocytes")
head(rM0)
DimPlot(rM0)
rM0 <- NormalizeData(rM0,normalization.method = "LogNormalize",scale.factor = 10000)
head(rM0)
rM0 <- FindVariableFeatures(rM0,selection.method = "vst",nfeatures = 2000)
top10 <- head(VariableFeatures(rM0),10)
plot3 <- VariableFeaturePlot(rM0)
plot4 <- LabelPoints(plot = plot3,points = top10,repel = T)
plot3 + plot4
rM0 <- ScaleData(rM0,features = rownames(rM0))
rM0 <- RunPCA(rM0,features = VariableFeatures(object = rM0))
print(rM0[["pca"]],dims = 1:5,nfeatures = 5)
VizDimLoadings(rM0,dims = 1:2,reduction = "pca")
DimPlot(rM0,reduction = "pca")
DimHeatmap(rM0,dims = 1,cells = 500,balanced = T)
DimHeatmap(rM0,dims = 1:15,cells = 500,balanced = T)
rM0 <- JackStraw(rM0,num.replicate = 100)
rM0 <- ScoreJackStraw(rM0,dims = 1:20)
JackStrawPlot(rM0,dims = 1:20)
ElbowPlot(rM0)
rM0 <- FindNeighbors(rM0,dims = 1:20)
rM0 <- FindClusters(rM0,resolution = 0.5)
rM0 <- RunUMAP(rM0,dims = 1:20)
head(rM0)
table(rM0@meta.data[["orig.ident"]])
DimPlot(rM0,reduction = "umap",pt.size = 0.8)
ggsave("monocyte plot.pdf",width = 10,height = 8)
M_markers0 <- FindAllMarkers(rM0,only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)
M_markers0 %>% group_by(cluster) %>% top_n(n = 2,wt = avg_log2FC)
M_top10m0 <- M_markers0 %>% group_by(cluster) %>% top_n(n = 10,wt = avg_log2FC)

VlnPlot(rM0,features = c("Aif1"),slot = "counts",log = T)
VlnPlot(rM0,features = M_top10m0[M_top10m0$cluster == 0,]$gene,slot = "counts",log = T)
FeaturePlot(rM0,features = c("Aif1"),
            cols = c("gray","red"),pt.size = 0.8) 

DotPlot(rM0,features = "Aif1")
ggsave("dotplot of aif1 in monocytes.pdf",width = 8,height = 6)
Idents(rM0)
deg <- FindMarkers(rM0,ident.1 = "1",ident.2 = "2",only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)
library(clusterProfiler)
library(org.Mm.eg.db)
gene <- rownames(deg)
gene.df <- bitr(gene, fromType = 'SYMBOL',
                toType = c('ENSEMBL','ENTREZID'),
                OrgDb = org.Mm.eg.db)
head(gene.df)
ee=enrichGO(gene = gene.df$ENTREZID,OrgDb = org.Mm.eg.db,pvalueCutoff = 0.05, ont="all", readable =T)
GO=as.data.frame(ee)
GO=GO[(GO$pvalue<0.05 & GO$qvalue<0.05),]
write.csv(GO,file="GO of diffgene in monocyte cluster 1 vs 2.csv")
showNum= 5
colorSel <- "qvalue"
pdf(file="GObarplot of diffgene in monocyte cluster 1 vs 2.pdf",width = 8,height =6)
bar=barplot(ee,drop = TRUE, showCategory =showNum,split="ONTOLOGY",color = colorSel)+ facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()

###KEGG enrichment analysis of DEGs
R.utils::setOption("clusterProfiler.download.method","auto")
kk <- enrichKEGG(gene = gene.df$ENTREZID,
                 organism = "mmu",
                 pvalueCutoff = 0.05)
KEGG=as.data.frame(kk)
KEGG
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(gene.df$SYMBOL[match(strsplit(x,"/")[[1]],as.character(gene.df$ENTREZID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<1 & KEGG$qvalue<1),]
write.csv(KEGG,file="KEGG of diffgene in monocyte cluster 1 vs 2.csv")
showNum=10

pdf(file="KEGGdotplot of in monocyte cluster 1 vs 2.pdf",width = 10,height = 8)
dotplot(kk,  orderBy = "GeneRatio", showCategory = showNum, color = colorSel)
dev.off()

library(SingleR)
library(celldex)
meta <- rM0@meta.data
immune_mus_gene <- celldex::ImmGenData()
exp_singleR<- GetAssayData(rM0, slot="data") 
m2.hpca <- SingleR(test = exp_singleR,ref = immune_mus_gene,labels = immune_mus_gene$label.fine)
head(m2.hpca)
rM0@meta.data$labels <- m2.hpca$labels
head(rM0)
DimPlot(rM0,reduction = "umap",pt.size = 0.5)
rM0@meta.data$cluseter <- Idents(rM0)
head(rM0)
head(r3@meta.data)
meta <- rM0@meta.data
DimPlot(rM0, group.by = c("labels"),reduction = "umap",split.by = "orig.ident",pt.size = 0.8)
