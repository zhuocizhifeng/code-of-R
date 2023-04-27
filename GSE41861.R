rm(list=ls())
library(GEOquery)
library(stringr)
library(data.table)
setwd("D:\\R\\asthma\\GSE41861")
r1 <- getGEO("GSE41861",getGPL = F)
r2 <- r1[[1]]
r3 <- pData(r2)
r4 <- as.data.frame(exprs(r2))
library(hgu133plus2.db)
gpl1 <- toTable(hgu133plus2SYMBOL)
ID <- rownames(r4)
r4$probe_id <- ID
eset1 <- merge(r4,gpl1,by="probe_id")
colnames(eset1)
eset2 <- aggregate(eset1[,2:139],by=list(eset1$symbol),FUN= mean)
rownames(eset2) <- eset2$Group.1
eset2 <- eset2[,-1]
boxplot(eset2)
write.csv(eset2,"Total exp.csv")
write.csv(r3,"clinical data.csv")
#upper vs control
group <- r3[r3$`tissue:ch1` == "Upper airway (Nasal)",]
nomalsample <- group[group$`disease:ch1` == "Control",]$geo_accession
assample <- group[group$`disease:ch1` == "Asthma",]$geo_accession
exp1 <- eset2[,nomalsample]
exp2 <- eset2[,assample]
eset3 <- cbind(exp1,exp2)
grouplist <- c(rep("CTL",17),rep("AA",40))
head(grouplist)
library(limma)
design <- model.matrix(~0+factor(grouplist))
colnames(design) <- levels(factor(grouplist))
contrast.matrix <- makeContrasts("AA-CTL",levels = design)
fit <- lmFit(eset3,design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
DEG <- topTable(fit2,number = Inf,adjust.method = "fdr")
head(DEG,20)
write.csv(DEG,"DEG of AA vs CTL in Nasal.csv")
diffgene <- DEG[with(DEG,(abs(logFC)>1 & adj.P.Val < 0.05)),]
head(diffgene)
write.csv(diffgene,file = "diffgene of AA vs CTL in Nasal.csv")

#lower vs control
group <- r3[r3$`tissue:ch1` == "Lower airway (bronchial)",]
nomalsample <- group[group$`disease:ch1` == "Control",]$geo_accession
assample <- group[group$`disease:ch1` == "Asthma",]$geo_accession
exp1 <- eset2[,nomalsample]
exp2 <- eset2[,assample]
eset3 <- cbind(exp1,exp2)
grouplist <- c(rep("CTL",30),rep("AA",51))
head(grouplist)
library(limma)
design <- model.matrix(~0+factor(grouplist))
colnames(design) <- levels(factor(grouplist))
contrast.matrix <- makeContrasts("AA-CTL",levels = design)
fit <- lmFit(eset3,design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
DEG <- topTable(fit2,number = Inf,adjust.method = "fdr")
head(DEG,20)
write.csv(DEG,"DEG of AA vs CTL in bronchial.csv")
diffgene <- DEG[with(DEG,(abs(logFC)>1 & adj.P.Val < 0.05)),]
head(diffgene)
write.csv(diffgene,file = "diffgene of AA vs CTL in bronchial.csv")

#Total asthma vs control
group <- r3
nomalsample <- group[group$`disease:ch1` == "Control",]$geo_accession
assample <- group[group$`disease:ch1` == "Asthma",]$geo_accession
exp1 <- eset2[,nomalsample]
exp2 <- eset2[,assample]
eset3 <- cbind(exp1,exp2)
grouplist <- c(rep("CTL",47),rep("AA",91))
head(grouplist)
library(limma)
design <- model.matrix(~0+factor(grouplist))
colnames(design) <- levels(factor(grouplist))
contrast.matrix <- makeContrasts("AA-CTL",levels = design)
fit <- lmFit(eset3,design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
DEG <- topTable(fit2,number = Inf,adjust.method = "fdr")
head(DEG,20)
write.csv(DEG,"DEG of AA vs CTL in Total.csv")
diffgene <- DEG[with(DEG,(abs(logFC)>1 & adj.P.Val < 0.05)),]
head(diffgene)
write.csv(diffgene,file = "diffgene of AA vs CTL in Total.csv")

#upper asthma vs lower asthma
group <- r3[r3$`tissue:ch1` == "Upper airway (Nasal)",]
uppersample <- group[group$`disease:ch1` == "Asthma",]$geo_accession
group <- r3[r3$`tissue:ch1` == "Lower airway (bronchial)",]
lowersample <- group[group$`disease:ch1` == "Asthma",]$geo_accession
exp1 <- eset2[,uppersample]
exp2 <- eset2[,lowersample]
eset3 <- cbind(exp1,exp2)
grouplist <- c(rep("UP",40),rep("LOWER",51))
head(grouplist)
library(limma)
design <- model.matrix(~0+factor(grouplist))
colnames(design) <- levels(factor(grouplist))
contrast.matrix <- makeContrasts("UP-LOWER",levels = design)
fit <- lmFit(eset3,design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
DEG <- topTable(fit2,number = Inf,adjust.method = "fdr")
head(DEG,20)
write.csv(DEG,"DEG of Asthma up vs lower.csv")
diffgene <- DEG[with(DEG,(abs(logFC)>1 & adj.P.Val < 0.05)),]
head(diffgene)
write.csv(diffgene,file = "diffgene of Asthma up vs lower.csv")
