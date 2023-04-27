library(data.table)
library(DESeq2)
setwd("D:\\R\\asthma\\GSE133278")
library(stringr)
wt1 <- fread("GSM3904690_PC-BM-1WT.1_sorted.counts.txt")
wt2 <- fread("GSM3904691_PC-BM-2WT.1_sorted.counts.txt")
wt3 <- fread("GSM3904692_PC-BM-3WT.1_sorted.counts.txt")
ko1 <- fread("GSM3904693_PC-BM-1KO.1_sorted.counts.txt")
ko2 <- fread("GSM3904694_PC-BM-2KO.1_sorted.counts.txt")
ko3 <- fread("GSM3904695_PC-BM-3KO.1_sorted.counts.txt")
eset1 <- data.frame(wt1 <- wt1$V2,
                    wt2 <- wt2$V2,
                    wt3 <- wt3$V2,
                    ko1 <- ko1$V2,
                    ko2 <- ko2$V2,
                    ko3 <- ko3$V2)
eset1$SYMBOL <- wt1$V1
library(org.Mm.eg.db)
gene <- as.data.frame(eset1$SYMBOL)
library(clusterProfiler)
gene1 <- bitr(gene$`eset1$SYMBOL`,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = org.Mm.eg.db)
head(gene1)
colnames(eset1)[7] <- "ENSEMBL"
eset2 <- merge(eset1,gene1,by = "ENSEMBL")
colnames(eset2)
eset2 <- aggregate(eset2[2:7],by=list(eset2$SYMBOL),FUN = mean)
rownames(eset2) <- eset2$Group.1
eset2 <- eset2[,-1]
colnames(eset2) <- c("GSM3904690","GSM3904691","GSM3904692",
                     "GSM3904693","GSM3904694","GSM3904695")
eset2 <- round(eset2)
grouplist <- c(rep("ctl",3),rep("ko",3))
condition <- factor(grouplist,levels = c("ctl","ko"))
colData <- data.frame(row.names = colnames(eset2),condition)
dds <- DESeqDataSetFromMatrix(countData = eset2, colData = colData, design = ~ condition)
dds <- dds[rowSums(counts(dds)) > 1,]
dds <- DESeq(dds)
deseq_normal_vst_expression <- assay(vst(dds,blind = T))
write.csv(deseq_normal_vst_expression,"gene exp after DESEQ2 normolized.csv")
res <- results(dds,contrast = c("condition","ko","ctl"))
resorder <- res[order(res$pvalue),]
DEG <- as.data.frame(resorder)
DEG <- na.omit(DEG)
logfc_cutoff <- 1
logfc_cutoff
DEG$change <- as.factor(ifelse(DEG$padj< 0.05 & abs(DEG$log2FoldChange)> logfc_cutoff,
                               ifelse(DEG$log2FoldChange > logfc_cutoff,"UP","DOWN"),"NOT"))
write.csv(DEG,"DEG of ctl vs aif1_ko.csv")
diffgene <- DEG[with(DEG,(abs(DEG$log2FoldChange)>1 & DEG$padj < 0.05)),]
head(diffgene)
write.csv(diffgene,"diffgene of ctl vs aif1ko.csv")

library(IOBR)
eset3 <- count2tpm(eset2,idType = "SYMBOL",org = "mmu")
eset3 <- log2(eset3+1)
eset3 <- as.data.frame(eset3)
write.csv(eset3,"gene exp of log2(TPM+1).csv")
eset4 <- as.data.frame(t(eset3))

library(homologene)
gene_red <- fread("D:\\R\\asthma\\GSE143303\\unique gene of red module among asthma group.csv")
gene_red <- gene_red$`gene symbol`
gene_red <- homologene(genes = gene_red,inTax = "9606",outTax = "10090")
genelist <- gene_red$`10090`
data1 <- eset4[,genelist[c(1:5,7:10)]]
data1$group <- data$group
gene_name <- colnames(data1)
for (i in 1:9) {
  gene1 <- gene_name[i]
  filename <- paste0(gene1,"_expression.pdf")
  p1 <- ggboxplot(data = data1,x="group",y= gene1,bxp.errorbar = T,
                  bxp.errorbar.width = 0.2,xlab = "Aif1 group",
                  ylab = "gene expression (log2 TPM+1)",ggtheme = theme_light(),
                  add = "jitter",fill = "group",palette = "lancet")
  p1 + stat_compare_means(method = "t.test")
  ggsave(filename = filename,width = 8,height = 6)
  
}

data2 <- eset4[,genelist[c(1:5,7,9,10)]]
dat <- data2
result <- data.frame(v1 = 1:(4 * nrow(dat)))  ## ?????ɺϲ????Ŀ??ܣ?��?кϲ?Ϊ1?У?????????????
for (i in 1:(ncol(dat)/4)) {                  ##??????????һ??????ѭ??
  temp1 <- c(dat[,i * 4 - 3], dat[, i * 4 - 2],
             dat[,i * 4 - 1], dat[, i * 4 - 0])   ##?ϲ?
  result <- cbind(result, temp1)              ##???ӵ????ݿ?
}
result
result <- result[,-1]
dat <- result
result <- data.frame(v1 = 1:(2 * nrow(dat)))  ## ?????ɺϲ????Ŀ??ܣ?��?кϲ?Ϊ1?У?????????????

for (i in 1:(ncol(dat)/2)) {                  ##??????????һ??????ѭ??
  temp1 <- c(dat[,i * 2 - 1], dat[, i * 2])   ##?ϲ?
  result <- cbind(result, temp1)              ##???ӵ????ݿ?
}
result
result <- result[,-1]
result <- as.data.frame(result)
genelist <- as.vector(colnames(data2))
op <- as.character()
for (i in genelist) {
  name <- i
  o <- rep(name,6)
  op <- append(op,o)
}
result$gene <- op
result$group <- c(rep(grouplist,8))
gp <- c(rep(grouplist,8))
combat_data <- result

p1 <- ggbarplot(data = combat_data,x= "gene",y = "result",fill = "group",
                palette = "lancet",position = position_dodge(),add = c("mean_se","point"),
                xlab = "red module gene",ylab = "gene expression (log2TPM+1)"
                )
p1 + stat_compare_means(aes(group = gp),label = "p.signif",method = "t.test")
ggsave(filename = "red module gene expression.pdf",width = 10,height = 8)


data <- as.data.frame(eset4$Il1b)
data$Ifnb1 <- eset4$Ifnb1
data$Il6 <- eset4$Il6
data$Il5 <- eset4$Il5
data$Il13 <- eset4$Il13
data$Tnf <- eset4$Tnf
data$IL10 <- eset4$Il10
data$Tgfb1 <- eset4$Tgfb1
colnames(data)[1] <- "Il1b"
dat <- data
result <- data.frame(v1 = 1:(4 * nrow(dat)))  ## ?????ɺϲ????Ŀ??ܣ?��?кϲ?Ϊ1?У?????????????
for (i in 1:(ncol(dat)/4)) {                  ##??????????һ??????ѭ??
  temp1 <- c(dat[,i * 4 - 3], dat[, i * 4 - 2],
             dat[,i * 4 - 1], dat[, i * 4 - 0])   ##?ϲ?
  result <- cbind(result, temp1)              ##???ӵ????ݿ?
}
result
result <- result[,-1]
dat <- result
result <- data.frame(v1 = 1:(2 * nrow(dat)))  ## ?????ɺϲ????Ŀ??ܣ?��?кϲ?Ϊ1?У?????????????

for (i in 1:(ncol(dat)/2)) {                  ##??????????һ??????ѭ??
  temp1 <- c(dat[,i * 2 - 1], dat[, i * 2])   ##?ϲ?
  result <- cbind(result, temp1)              ##???ӵ????ݿ?
}
result
result <- result[,-1]
result <- as.data.frame(result)
genelist <- as.vector(colnames(data))
op <- as.character()
for (i in genelist) {
  name <- i
  o <- rep(name,6)
  op <- append(op,o)
}
result$gene <- op
result$group <- c(rep(grouplist,8))
gp <- c(rep(grouplist,8))
combat_data <- result
library(ggpubr)
p1 <- ggbarplot(data = combat_data,x= "gene",y = "result",fill = "group",
                palette = "lancet",position = position_dodge(),add = c("mean_se","point"),
                xlab = "Macrophagy related cytokine",ylab = "gene expression (log2TPM+1)"
)
p1 + stat_compare_means(aes(group = gp),label = "p.signif",method = "t.test")
ggsave(filename = "macrophagy cytokine.pdf",width = 10,height = 8)
eset4$Ly6c1
chemogene <- c("Ccl2","Ccl8","Ccl17","Ccl22","Cxcl1","Cxcl2","Cxcl9","Cxcl10","Cxcl11","Ccl5")

data_c <- eset4[,chemogene]
dat <- data_c
result <- data.frame(v1 = 1:(5 * nrow(dat)))  ## ?????ɺϲ????Ŀ??ܣ?��?кϲ?Ϊ1?У?????????????
for (i in 1:(ncol(dat)/5)) {                  ##??????????һ??????ѭ??
  temp1 <- c(dat[,i * 5 - 4], dat[, i * 5 - 3],
             dat[,i * 5 - 2], dat[, i * 5 - 1],dat[, i * 5 - 0])   ##?ϲ?
  result <- cbind(result, temp1)              ##???ӵ????ݿ?
}
result
result <- result[,-1]
dat <- result
result <- data.frame(v1 = 1:(2 * nrow(dat)))  ## ?????ɺϲ????Ŀ??ܣ?��?кϲ?Ϊ1?У?????????????
for (i in 1:(ncol(dat)/2)) {                  ##??????????һ??????ѭ??
  temp1 <- c(dat[,i * 2 - 1], dat[, i * 2])   ##?ϲ?
  result <- cbind(result, temp1)              ##???ӵ????ݿ?
}
result
result <- result[,-1]
result <- as.data.frame(result)
genelist <- as.vector(colnames(data_c))
op <- as.character()
for (i in genelist) {
  name <- i
  o <- rep(name,6)
  op <- append(op,o)
}
result$gene <- op
result$group <- c(rep(grouplist,10))
gp <- c(rep(grouplist,10))
combat_data <- result
library(ggpubr)
p1 <- ggbarplot(data = combat_data,x= "gene",y = "result",fill = "group",
                palette = "lancet",position = position_dodge(),add = c("mean_se","point"),
                xlab = "Macrophagy related chemokines",ylab = "gene expression (log2TPM+1)"
)
p1 + stat_compare_means(aes(group = gp),label = "p.signif",method = "t.test")
ggsave(filename = "macrophagy chemokine.pdf",width = 10,height = 8)

M1_gene <- c("Il1a","Il1b","Il6","Nos2","Tlr2","Tlr4","Cd80","Cd86","Tnf","Ly6c1")
M2_gene <- c("Cd163","Pparg","Arg1","Csf1r","Mrc1","Il10")
M_gene <- c("Il1a","Il1b","Il6","Nos2","Tlr2","Tlr4","Cd80","Cd86","Tnf","Ly6c1",
            "Cd163","Pparg","Arg1","Csf1r","Mrc1","Il10")

data_c <- eset4[,M_gene]
dat <- data_c
result <- data.frame(v1 = 1:(8 * nrow(dat)))  
for (i in 1:(ncol(dat)/8)) {                  
  temp1 <- c(dat[,i * 8 - 7], dat[, i * 8 - 6],dat[,i * 8 - 5], dat[,i * 8 - 4],dat[, i * 8 - 3],
             dat[,i * 8 - 2], dat[, i * 8 - 1],dat[, i * 8 - 0])  
  result <- cbind(result, temp1)             
}
result
result <- result[,-1]
dat <- result
result <- data.frame(v1 = 1:(2 * nrow(dat)))  ## ?????ɺϲ????Ŀ??ܣ?��?кϲ?Ϊ1?У?????????????
for (i in 1:(ncol(dat)/2)) {                  ##??????????һ??????ѭ??
  temp1 <- c(dat[,i * 2 - 1], dat[, i * 2])   ##?ϲ?
  result <- cbind(result, temp1)              ##???ӵ????ݿ?
}
result
result <- result[,-1]
result <- as.data.frame(result)
genelist <- as.vector(colnames(data_c))
op <- as.character()
for (i in genelist) {
  name <- i
  o <- rep(name,6)
  op <- append(op,o)
}
result$gene <- op
result$group <- c(rep(grouplist,16))
gp <- c(rep(grouplist,16))
combat_data <- result
library(ggpubr)
p1 <- ggbarplot(data = combat_data,x= "gene",y = "result",fill = "group",
                palette = "lancet",position = position_dodge(),add = c("mean_se","point"),
                xlab = "Macrophagy Marker",ylab = "gene expression (log2TPM+1)"
)
p1 + stat_compare_means(aes(group = gp),label = "p.signif",method = "t.test")
ggsave(filename = "macrophagy M1 or M2.pdf",width = 10,height = 8)

















library(ggplot2)
p1 <- ggboxplot(data = data1,x="group",y= "Ly96",bxp.errorbar = T,
                bxp.errorbar.width = 0.2,xlab = "Aif1 group",
                ylab = "gene expression (log2 TPM+1)",ggtheme = theme_light(),
                add = "jitter",fill = "group",palette = "lancet")
p1 + stat_compare_means()
ggsave(filename = "Il4 expression.pdf",width = 8,height = 6)

head(DEG)
head(diffgene)
#GO for DEGs
library(clusterProfiler)
library(org.Mm.eg.db)
colorSel="qvalue"
gene <- rownames(diffgene)
#entrezIDs <- mget(gene, org.Hs.egSYMBOL2EG, ifnotfound=NA)
gene.df <- bitr(gene, fromType = 'SYMBOL',
                toType = c('ENSEMBL','ENTREZID'),
                OrgDb = org.Mm.eg.db)
head(gene.df)
ee=enrichGO(gene = gene.df$ENTREZID,OrgDb = org.Mm.eg.db,ont="all", readable =T,
            qvalueCutoff = 0.05)
GO=as.data.frame(ee)
write.csv(GO,file="GO_DEGs.csv")
library(ggplot2)
pdf(file="DEGs_GObarplot.pdf",width = 10,height =8)
bar=barplot(ee, drop = TRUE, showCategory = 5,split="ONTOLOGY",color = colorSel)+ facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()
#KEGG
library(DO.db)
R.utils::setOption("clusterProfiler.download.method","auto")
kk <- enrichKEGG(gene = gene.df$ENTREZID,
                 organism = "mmu",
                 pvalueCutoff = 0.05)
KEGG=as.data.frame(kk)
KEGG
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(gene.df$SYMBOL[match(strsplit(x,"/")[[1]],as.character(gene.df$ENTREZID))],collapse="/")))

write.csv(KEGG,file="KEGG of deg.csv")

pdf(file="DEGs_KEGGbarplot.pdf",width = 10,height = 8)
barplot(kk, drop = TRUE, showCategory = 10, color = colorSel)
dev.off()


DEG1 <- DEG[order(DEG$log2FoldChange,decreasing = T),]
genelist <- DEG1$log2FoldChange
names(genelist) <- rownames(DEG1)
genelist1 <- as.data.frame(genelist)
genelist1$symbol <- rownames(genelist1)
gene1 <- bitr(genelist1$symbol,OrgDb = org.Mm.eg.db,fromType = "SYMBOL",toType = "ENTREZID")
colnames(gene1)[1] <- "symbol"
genelist2 <- merge(genelist1,gene1,by= "symbol")
genelist3 <- genelist2[order(genelist2$genelist,decreasing = T),]
genelist4 <- genelist3$genelist
names(genelist4) <- genelist3$ENTREZID
library(enrichplot)
gse.GO <- gseGO(
  genelist4,
  ont = "ALL",  
  OrgDb = org.Mm.eg.db, 
  keyType = "ENTREZID",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)
go <- as.data.frame(gse.GO)
go$core_enrichment <- as.character(sapply(go$core_enrichment,function(x)paste(genelist2$symbol[match(strsplit(x,"/")[[1]],as.character(genelist2$ENTREZID))],collapse="/")))
write.csv(go,"GSEA GO of aif1 ko.csv")
pdf("GSEA GO dotplot top10 of aif1 ko.pdf",width=10,height = 10)
dotplot(gse.GO)
dev.off()
pdf("GSEA go GSEAplot top10 of aif1 ko.pdf",width = 10,height = 10)
p <- gseaplot2(gse.GO,1:10)
print(p)
dev.off()
pdf("GSEA GO ridgeplot top10 of aif1 ko.pdf",width=10,height = 10)
ridgeplot(gse.GO,showCategory = 10)
dev.off()
gse.kegg <- gseKEGG(genelist4,organism = "mmu",
                    pvalueCutoff = 0.05)
kegg <- as.data.frame(gse.kegg)
kegg$core_enrichment <- as.character(sapply(kegg$core_enrichment,function(x)paste(genelist2$symbol[match(strsplit(x,"/")[[1]],as.character(genelist2$ENTREZID))],collapse="/")))
write.csv(kegg,"GSEA kegg of 10ug vs 100ug.csv")
pdf("GSEA KEGG dotplot top10 of aif1 ko.pdf",width=10,height = 10)
dotplot(gse.kegg)
dev.off()
pdf("GSEA KEGG GSEAplot of aif1 ko.pdf",width = 10,height = 10)
p <- gseaplot2(gse.kegg,1:5)
print(p)
dev.off()
pdf("GSEA kegg ridgeplot aif1 ko.pdf",width=10,height = 10)
ridgeplot(gse.kegg,showCategory = 10)
dev.off()

data <- deseq_normal_vst_expression
library(GSVA)
library(GSEABase)
library(limma)
library(stringr)
geneset <- getGmt("mh.all.v2022.1.Mm.symbols.gmt")
data1 <- as.matrix(data)
es <- gsva(data1,geneset,verbose=T,mx.diff=T)
rownames(es) <- str_replace(rownames(es),"HALLMARK_","")
exp4 <- as.data.frame(t(data))
design <- model.matrix(~0+factor(grouplist))
colnames(design) <- levels(factor(grouplist))
grouplist
contrast.matrix <- makeContrasts("ko-ctl",levels = design)
fit <- lmFit(es,design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
DEG2 <- topTable(fit2,number = Inf,adjust.method = "fdr")
head(DEG2,20)
write.csv(DEG2,"DEpathway of aif1 ko.csv")
diffgene2 <- DEG2[with(DEG2,(abs(t)>2 & adj.P.Val < 0.05)),]
head(diffgene2)
write.csv(diffgene2,file = "diff pathway of aif1 ko.csv")

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
ggsave("gsva_pathway_bar_aif1 ko.pdf",p,width = 14,height  = 10)

geneset <- read.csv("genesets.csv")
set1 <- geneset[1]
set1 <- as.data.frame(set1)
gene <- set1$myogenesis
data_t <- as.data.frame(t(data))
id <- as.data.frame(colnames(data_t))
gene1 <- as.data.frame(gene) 
colnames(id) <- "gene"
gene2 <- merge(gene1,id,by = "gene")
choose_gene = gene2$gene
write.csv(choose_gene,"myogenesis gene.csv")
choose_matrix = data[choose_gene,]
choose_matrix = t(scale(t(choose_matrix)))
annotation_col=data.frame(group=grouplist)
rownames(annotation_col)=colnames(data)
pdf("heatmap of myogenesis gene.pdf",height = 18,width = 16)
p = pheatmap(choose_matrix,color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
             scale = "row",annotation_col=annotation_col,show_rownames = T,show_colnames = F,
             main = "TOP50 heatmap ",cluster_rows = T,cluster_cols = F)
p
dev.off()

data_plot <- as.data.frame(eset4$Tnf)
head(data_plot)
colnames(data_plot) <- "Tnf"
data_plot$group <- grouplist
data_plot$Sting <- eset4$Sting1
data_plot$FCN1 <- data_nea$FCN1
data_plot$IL6<- eset4$Il6
data_plot$TNFR1B <- data_nea$TNFRSF1B
library(ggpubr)
p1 <- ggbarplot(data = data_plot,x= "group",y = "IL6",fill = "group",
                palette = "lancet",position = position_dodge(),add = c("mean_se","point"),
                xlab = "Aif1 status",ylab = "asc gene expression (log2TPM+1)",
                size = 1,ggtheme = theme_light())
p1 + stat_compare_means(method = "t.test")
ggsave("ASC expression.pdf",width=8,height = 6)

