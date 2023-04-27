library(Seurat)
library(data.table)
library(ggplot2)
setwd("D:\\R\\asthma\\GSE155391")
r1 <- fread("GSE155391_umi-count-matrix.csv")
ID <- r1$genes
r1 <- r1[,-1]
r1 <- r1[,c(1:4912)]
tail(colnames(r1))
rownames(r1) <- ID
rownames(r1)
colnames(r1)
r3 <- fread("GSE155391_cell-level-metadata.csv")
table(r3$treatment)
r3 <- r3[c(1:4912),]
r2 <- CreateSeuratObject(counts = r1)
rownames(r2)
r2 <- AddMetaData(object = r2,metadata = r3)
r2[["percent.mt"]] <- PercentageFeatureSet(r2,pattern = "^mt-")
head(r2)
dim(r2)
r2$treatment <- c(rep("control",3918),rep("HDM",994))
table(r2$treatment)
colnames(r2)
library(tidyverse)
VlnPlot(r2,features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3)
plot1 <- FeatureScatter(r2,feature1 = "nCount_RNA",feature2 = "percent.mt")
plot2 <- FeatureScatter(r2,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1,plot2))
r2 <- subset(r2,subset = nFeature_RNA > 250 & nFeature_RNA < 4000 & percent.mt < 10)
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
r3 <- FindNeighbors(r3,dims = 1:18)
r3 <- FindClusters(r3,resolution = 0.8)
r3 <- RunUMAP(r3,dims = 1:18)
r3 <- RunTSNE(r3,dims = 1:18)
DimPlot(r3,reduction = "umap",pt.size = 0.8,split.by = "treatment")
DimPlot(r3,reduction = "umap",pt.size = 0.8)
markers <- FindAllMarkers(r3,only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)
markers %>% group_by(cluster) %>% top_n(n = 2,wt = avg_log2FC)
VlnPlot(r3,features = c("Il1b"),slot = "counts",log = T)
FeaturePlot(r3,features = c("LYZ","RUNX2","CD3D","CTSK","RGS5","TOP2A","EGFL7","ACAN"),
            cols = c("gray","red")) 
ggsave("Figure 1E.pdf",width = 12,height = 12)
top10m <- markers %>% group_by(cluster) %>% top_n(n = 10,wt = avg_log2FC)
DoHeatmap(r3,features = top10m$gene) + NoLegend()
top10m[top10m$cluster == 0,]$gene
VlnPlot(r3,features = top10m[top10m$cluster == 11,]$gene,slot = "counts",log = T)
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
r4$cell <- Idents(r3)
r4$cell.group <- paste(r3$cell,r3$treatment,sep = "_")
unique(r3$cell.group)
Idents(r4) <- "cell.group"
r3$cell.group
cellfordeg <- levels(r3$cell)
for (i in c(1,2,4,5,6)) {
  CELLDEG <- FindMarkers(r4,ident.1 = paste0(cellfordeg[i],"_control"),
                         ident.2 = paste0(cellfordeg[i],"_HDM"),
                         verbose = F,
                         test.use = "wilcox")
  write.csv(CELLDEG,paste0(cellfordeg[i],".csv"))
  
}
r4$treatment
Idents(r4) <- "treatment"
VlnPlot(r4,features = c("Hmox1","Lyz2","Mt1","Mt2"))
FeaturePlot(r4,features = c("Traf4"))
DotPlot(r3,features = c("Hmox1","Lyz2","Mt1","Mt2"),group.by = "treatment")
