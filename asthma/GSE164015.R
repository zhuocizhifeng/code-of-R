setwd("D:\\R\\asthma\\GSE164015")
library(Seurat)
library(dplyr)
library(multtest)
library(mindr)
library(tidyverse)
library(GEOquery)
library(data.table)
r1 <- getGEO("GSE164015",getGPL = F)
r2 <- r1[[1]]
r3 <- pData(r2)
c1 <- Read10X(data.dir = "D:\\R\\asthma\\GSE164015\\ACE043DIL_Epi_Raw") 
c1 <- CreateSeuratObject(counts = c1, project = "Control",min.cells = 3,min.features = 200)
head(c1)
c1[["percent.mt"]] <- PercentageFeatureSet(c1,pattern = "^MT-")
c2 <- Read10X(data.dir = "D:\\R\\asthma\\GSE164015\\ACE049DIL_Epi_Raw") 
c2 <- CreateSeuratObject(counts = c2, project = "Control",min.cells = 3,min.features = 200)
head(c2)
c2[["percent.mt"]] <- PercentageFeatureSet(c2,pattern = "^MT-")
c3 <- Read10X(data.dir = "D:\\R\\asthma\\GSE164015\\ACE057DIL_Epi_Raw") 
c3 <- CreateSeuratObject(counts = c3, project = "Control",min.cells = 3,min.features = 200)
head(c3)
c3[["percent.mt"]] <- PercentageFeatureSet(c3,pattern = "^MT-")
t1 <- Read10X(data.dir = "D:\\R\\asthma\\GSE164015\\ACE043AC_Epi_Raw") 
t1 <- CreateSeuratObject(counts = t1, project = "HDM",min.cells = 3,min.features = 200)
head(t1)
t1[["percent.mt"]] <- PercentageFeatureSet(t1,pattern = "^MT-")
t2 <- Read10X(data.dir = "D:\\R\\asthma\\GSE164015\\ACE049AC_Epi_Raw") 
t2 <- CreateSeuratObject(counts = t2, project = "HDM",min.cells = 3,min.features = 200)
head(t2)
t2[["percent.mt"]] <- PercentageFeatureSet(t2,pattern = "^MT-")
t3 <- Read10X(data.dir = "D:\\R\\asthma\\GSE164015\\ACE057AC_Epi_Raw") 
t3 <- CreateSeuratObject(counts = t3, project = "HDM",min.cells = 3,min.features = 200)
head(t3)
t3[["percent.mt"]] <- PercentageFeatureSet(t3,pattern = "^MT-")
list1 <- list(c2,c3,t1,t2,t3)
m1 <- merge(x=c1,y=list1)
head(m1)
VlnPlot(m1,features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3)
m2 <- subset(m1,subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
VlnPlot(m2,features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3)
m2 <- NormalizeData(m2,normalization.method = "LogNormalize",scale.factor = 10000)
m2 <- FindVariableFeatures(m2,selection.method = "vst",nfeatures = 2000)
top10 <- head(VariableFeatures(m2),10)
plot3 <- VariableFeaturePlot(m2)
plot4 <- LabelPoints(plot = plot3,points = top10,repel = T)
plot3 + plot4
m2 <- ScaleData(m2,features = rownames(m2))
m2 <- RunPCA(m2,features = VariableFeatures(object = m2))
print(m2[["pca"]],dims = 1:5,nfeatures = 5)
VizDimLoadings(m2,dims = 1:2,reduction = "pca")
DimPlot(m2,reduction = "pca")
DimHeatmap(m2,dims = 1,cells = 500,balanced = T)
DimHeatmap(m2,dims = 1:15,cells = 500,balanced = T)
m2 <- JackStraw(m2,num.replicate = 100)
m2 <- ScoreJackStraw(m2,dims = 1:20)
JackStrawPlot(m2,dims = 1:20)
ElbowPlot(m2)
m2 <- FindNeighbors(m2,dims = 1:20)
m2 <- FindClusters(m2,resolution = 0.1)
m2 <- RunUMAP(m2,dims = 1:20)
head(m2)
dim(m2)
table(m2$orig.ident)
DimPlot(m2,reduction = "umap",pt.size = 0.5,split.by = "orig.ident")
markers_df <- FindMarkers(object = m2, ident.1 = 0, min.pct = 0.25)
markers_genes =  rownames(head(x = markers_df, n = 5))
markers <- FindAllMarkers(m2,only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)
top10m <- markers %>% group_by(cluster) %>% top_n(n = 10,wt = avg_log2FC)

genes <- markers[markers$cluster == 0,]$gene
VlnPlot(m2,features =top10m[top10m$cluster == 0,]$gene,slot = "counts",log = T)
library(COSG)
marker_cosg <- cosg(
  m2,
  groups='all',
  assay='RNA',
  slot='data',
  mu=1,
  n_genes_user=100)
g1 <- head(marker_cosg$names$`0`,10)
th= theme(axis.text.x = element_text(angle = 90))
g2 <- head(top10m$gene,10)
DotPlot(m2, features = g2,assay='RNA'  ) + th +NoLegend()
FeaturePlot(r3,features = c("LYZ","RUNX2","CD3D","CTSK","RGS5","TOP2A","EGFL7","ACAN"),
            cols = c("gray","red"))
#chang shi singleR
library(SingleR)
library(celldex)
hpca.se <- HumanPrimaryCellAtlasData()
exp_singleR<- GetAssayData(m2, slot="data") 
m2.hpca <- SingleR(test = exp_singleR,ref = hpca.se,labels = hpca.se$label.main)
head(m2.hpca)
m2@meta.data$labels <- m2.hpca$labels
head(m2)
DimPlot(m2,reduction = "umap",pt.size = 0.5)
m2@meta.data$cluseter <- Idents(m2)
head(m2)
head(m2@meta.data)
meta <- m2@meta.data
DimPlot(m2, group.by = c("labels"),reduction = "umap",split.by = "orig.ident")
#cell marker database
VlnPlot(m2,features =top10m[top10m$cluster == 1,]$gene,slot = "counts",log = T)
VlnPlot(m2,features = "CCL23",slot = "counts",log = T)
#0.Epithelial cell SCGB3A1
#1.Ciliated Cell CASP C20orf85 TUBB4B
#2.Ciliated Cell
#3.Secretory cell
#4.Basal cell
#5.Basal cell
#6.
#7.Ciliated Cell TMEM190
new.cluster.ids <- c("Epithelial_cells",
                     "Epithelial_cells",
                     "T cells",
                     )
names(new.cluster.ids) <- levels(r3)
r3 <- RenameIdents(r3,new.cluster.ids)
DimPlot(r3,reduction = "umap",label = T,pt.size = 0.5)
