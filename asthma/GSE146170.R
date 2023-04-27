library(Seurat)
library(data.table)
library(ggplot2)
setwd("D:\\R\\asthma\\GSE146170")
r1 <- fread("GSE146170_allergen_TH2_umi.txt")
ID <- r1$gene_name
r1 <- r1[,-1]
tail(colnames(r1))
rownames(r1) <- ID
rownames(r1)
colnames(r1)
meta <- fread("allergen_TH2_annotation.txt")
r2 <- CreateSeuratObject(counts = r1)
rownames(r2)
r2$percent.mt <- meta$percent.mt
r2$group <- meta$diseasegroup
head(r2)
dim(r2)
colnames(r2)
library(tidyverse)
VlnPlot(r2,features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3)
plot1 <- FeatureScatter(r2,feature1 = "nCount_RNA",feature2 = "percent.mt")
plot2 <- FeatureScatter(r2,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1,plot2))
r2 <- subset(r2,subset = nFeature_RNA > 250 & nFeature_RNA < 2000 & percent.mt < 5)
dim(r2)
r3 <- NormalizeData(r2,normalization.method = "LogNormalize",scale.factor = 10000)
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
r3 <- FindNeighbors(r3,dims = 1:13)
r3 <- FindClusters(r3,resolution = 0.5)
r3 <- RunUMAP(r3,dims = 1:13)
r3 <- RunTSNE(r3,dims = 1:13)
DimPlot(r3,reduction = "umap",pt.size = 0.8,split.by = "group")
DimPlot(r3,reduction = "umap",pt.size = 1)
markers <- FindAllMarkers(r3,only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)
markers %>% group_by(cluster) %>% top_n(n = 2,wt = avg_log2FC)
VlnPlot(r3,features = c("TXN"),slot = "counts",log = T)
FeaturePlot(r3,features = c("MT1X","MT2A","SOD2"),
            cols = c("gray","red")) 
ggsave("Figure 1E.pdf",width = 12,height = 12)
top10m <- markers %>% group_by(cluster) %>% top_n(n = 10,wt = avg_log2FC)
DoHeatmap(r3,features = top10m$gene) + NoLegend()
top10m[top10m$cluster == 0,]$gene
VlnPlot(r3,features = top10m[top10m$cluster == 4,]$gene,slot = "counts",log = T)
#cluster4 : cd3d,nkg7
#cluster8 : B cell Ly6d,cd79a
#cluster9 : neutrophil S100a9,MMP9 
#cluster10 : DC,fscn1,Nr4a3
#cluster11 : Ciliated cell Nupr1 Phlda1
ggsave("Figure 1D.pdf",width = 12,height = 12)
new.cluster.ids <- c("Alveolar macrophages",
                     "Alveolar macrophages",
                     "Alveolar macrophages",
                     "Alveolar macrophages",
                     "T cells",
                     "Alveolar macrophages",
                     "Alveolar macrophages",
                     "Alveolar macrophages",
                     "B cells",
                     "Neutrophils",
                     "Dendritic cells",
                     "Fibroblast")
names(new.cluster.ids) <- levels(r3)
r3 <- RenameIdents(r3,new.cluster.ids)
DimPlot(r3,reduction = "umap",label = T,pt.size = 0.5,split.by = "treatment",
        label.size = 2.5)
DimPlot(r3,reduction = "umap",pt.size = 0.8)
Idents(r3)
r4 <- r3
Idents(r4) <- r4$group
r4 <- subset(r4,idents = "AS_AL")
head(r4)
VlnPlot(r4,features = "MT2A")
head(r4$seurat_clusters)
r4 <- ScaleData(r4,features = rownames(r4))
r4 <- RunPCA(r4,features = VariableFeatures(object = r4))
r4 <- JackStraw(r4,num.replicate = 100)
r4 <- ScoreJackStraw(r4,dims = 1:20)
JackStrawPlot(r4,dims = 1:20)
ElbowPlot(r4)
r4 <- FindNeighbors(r4,dims = 1:13)
r4 <- FindClusters(r4,resolution = 0.5)
r4 <- RunUMAP(r4,dims = 1:13)
r4 <- RunTSNE(r4,dims = 1:13)
DimPlot(r4,reduction = "umap",pt.size = 0.8,split.by = "group")
DimPlot(r4,reduction = "tsne",pt.size = 1)
markers <- FindAllMarkers(r4,only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)
markers %>% group_by(cluster) %>% top_n(n = 2,wt = avg_log2FC)
top10m4 <- markers %>% group_by(cluster) %>% top_n(n = 10,wt = avg_log2FC)
DoHeatmap(r4,features = top10m4$gene) 
top10m4[top10m4$cluster == 1,]$gene
gene.name <- markers$gene[markers$cluster == "1"][1:20]
VlnPlot(r4,features = top10m4[top10m4$cluster == 3,]$gene,slot = "counts",log = T)
VlnPlot(r4,features = gene.name,slot = "counts",log = T)
VlnPlot(r4,features = c("IL13"),slot = "counts",log = T)
new.cluster.ids <- c("WDR43-","NOP14+","HILPDA+","SLC7A5")
names(new.cluster.ids) <- levels(r4)
r4 <- RenameIdents(r4,new.cluster.ids)
DimPlot(r4,reduction = "umap",label = T,pt.size = 1.5,
        label.size = 2.5,label.color = "black")
DimPlot(r4,reduction = "umap",pt.size = 0.8)