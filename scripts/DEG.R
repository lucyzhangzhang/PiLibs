library(DESeq2) #DEG analysis
library(dplyr) #data manipulation

wd <- "~/R/Eutrema/Pi/"
setwd(wd)

condition <- factor(rep(c("L0", "L2.5", "R0", "R2.5"), each = 2))
sample <- c("YLP0-1",
            "YLP0-2",
            "YLP2.5-1",
            "YLP2.5-2",
            "YRP0-1",
            "YRP0-2",
            "YRP2.5-1",
            "YRP2.5-2")

metadata <- data.frame(sampleNO = sample,
                       condition = condition,
                       libraryName = sample)

metadata <- mutate(metadata, countFile = paste0(metadata$libraryName, "/", metadata$libraryName,
                                                "QC.geneCounts.formatted.for.DESeq.txt"))

sampleTable <- data.frame(sampleName = metadata$libraryName,
                          fileName = metadata$countFile,
                          condition = metadata$condition,
                          sampleNO = metadata$sampleNO)


DESeq2Table <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                          directory = "QoRTs",
                                          design = ~ condition)

metadata$condition <- relevel(metadata$condition, ref = "L2.5")

keep <- rowSums(counts(DESeq2Table)) >= 8
DESeq2Table <- DESeq2Table[keep,]

dds <- DESeq(DESeq2Table)
resultsNames(dds)

ddsL2.5 <- dds
ddsL2.5$condition <- relevel(ddsL2.5$condition, ref = "L2.5")
ddsL2.5 <- nbinomWaldTest(ddsL2.5)
resultsNames(ddsL2.5)

L2.5_L0 <- results(ddsL2.5, contrast = c("condition", "L0", "L2.5"))
L2.5_R0 <- results(ddsL2.5, contrast = c("condition", "R0", "L2.5"))
L2.5_R2.5 <- results(ddsL2.5, contrast = c("condition", "R2.5", "L2.5"))

L2L0.up <- L2.5_L0[which(L2.5_L0[,"padj"] < 0.05 & L2.5_L0[,"log2FoldChange"] > 0),]
L2L0.dn <- L2.5_L0[which(L2.5_L0[,"padj"] < 0.05 & L2.5_L0[,"log2FoldChange"] < 0),]
L2R0.up <- L2.5_R0[which(L2.5_R0[,"padj"] < 0.05 & L2.5_R0[,"log2FoldChange"] > 0),]
L2R0.dn <- L2.5_R0[which(L2.5_R0[,"padj"] < 0.05 & L2.5_R0[,"log2FoldChange"] < 0),]
L2R2.up <- L2.5_R2.5[which(L2.5_R2.5[,"padj"] < 0.05 & L2.5_R2.5[,"log2FoldChange"] > 0),]
L2R2.dn <- L2.5_R2.5[which(L2.5_R2.5[,"padj"] < 0.05 & L2.5_R2.5[,"log2FoldChange"] < 0),]

nrow(L2L0.up)
nrow(L2L0.dn)
nrow(L2R0.up)
nrow(L2R0.dn)
nrow(L2R2.up)
nrow(L2R2.dn)
 write.table(L2L0.up, file = "data/L2L0.up", sep = "\t", row.names = T, col.names = F)
 write.table(L2L0.dn, file = "data/L2L0.dn", sep = "\t", row.names = T, col.names = F)

# redos the nbinomWalTest with the new relevelled reference
ddsR2.5 <- dds
ddsR2.5$condition <- relevel(ddsR2.5$condition, ref = "R2.5")
ddsR2.5 <- nbinomWaldTest(ddsR2.5)
resultsNames(ddsR2.5)

R2.5_L0 <- results(ddsR2.5, contrast = c("condition", "L0", "R2.5"))
R2.5_L2.5 <- results(ddsR2.5, contrast = c("condition", "L2.5", "R2.5"))
R2.5_R0 <- results(ddsR2.5, contrast = c("condition", "R0", "R2.5"))

R2L0.up <- R2.5_L0[which(R2.5_L0[,"padj"] < 0.05 & R2.5_L0[,"log2FoldChange"] > 0),]
R2L0.dn <- R2.5_L0[which(R2.5_L0[,"padj"] < 0.05 & R2.5_L0[,"log2FoldChange"] < 0),]
R2L2.up <- R2.5_L2.5[which(R2.5_L2.5[,"padj"] < 0.05 & R2.5_L2.5[,"log2FoldChange"] > 0),]
R2L2.dn <- R2.5_L2.5[which(R2.5_L2.5[,"padj"] < 0.05 & R2.5_L2.5[,"log2FoldChange"] < 0),]
R2R0.up <- R2.5_R0[which(R2.5_R0[,"padj"] < 0.05 & R2.5_R0[,"log2FoldChange"] > 0),]
R2R0.dn <- R2.5_R0[which(R2.5_R0[,"padj"] < 0.05 & R2.5_R0[,"log2FoldChange"] < 0),]

nrow(R2L0.up)
nrow(R2L0.dn)
nrow(R2L2.up)
nrow(R2L2.dn)
nrow(R2R0.up)
nrow(R2R0.dn)

 write.table(R2R0.up, file = "data/R2R0.up", sep = "\t", row.names = T, col.names = F)
 write.table(R2R0.dn, file = "data/R2R0.dn", sep = "\t", row.names = T, col.names = F)
L0_R2.5 <- results(dds, contrast = c("condition", "R2.5", "L0"))

L2R2.up <- L2.5_R2.5[which(L2.5_R2.5[,"padj"] < 0.05 & L2.5_R2.5[,"log2FoldChange"] > 0),]
L2R2.dn <- L2.5_R2.5[which(L2.5_R2.5[,"padj"] < 0.05 & L2.5_R2.5[,"log2FoldChange"] < 0),]

write.table(colnames(R2R0.up), file = "data/column.names", quote = F, row.names = F, sep = "\t")
nrow(L2R2.up)
nrow(L2R2.dn)


# Venn diagrams
library(GenomicFeatures)               #GFF functions makeTxdbFromGff
library(venn)

# get median length of libraries
# length correction
GTF <- makeTxDbFromGFF("~/R/Eutrema/PS/drought.gtf", format = "gtf")

GTF <- exonsBy(GTF, by = "gene")
head(GTF)
lengths <- sum(width(reduce(GTF)))
lengths <- as.data.frame(lengths)
inter2 <- intersect(rownames(lengths), rownames(dds))

length(inter2)
lens <- lengths[row.names(lengths) %in% inter2,]
lenMedian <- median(lens)

mcols(dds)$basepairs <- lenMedian
fpkm <- fpkm(dds, robust = T)
head(fpkm)
fpkm[row.names(fpkm) %in% "nXLOC_008023", ]

fpkm.ave <- data.frame(YLP0=rowMeans(fpkm[,c(1,2)]), YLP2.5=rowMeans(fpkm[,c(3,4)]),
                       YRP0=rowMeans(fpkm[,c(5,6)]), YRP2.5=rowMeans(fpkm[,c(7,8)]))

# perform set operation to see what genes are in each set
# e.g 
# setdiff(intersect(YRP2.5, YLP2.5), union(YRP0, YLP0))
# [1] "Thhalv10020696m.g"
YLP0 <- rownames(fpkm.ave[which(fpkm.ave$YLP0 != 0),])
YLP2.5 <- rownames(fpkm.ave[which(fpkm.ave$YLP2.5 != 0),])
YRP0 <- rownames(fpkm.ave[which(fpkm.ave$YRP0 != 0),])
YRP2.5 <- rownames(fpkm.ave[which(fpkm.ave$YRP2.5 != 0),])


PiSet <- list(YLP0=rownames(fpkm.ave[which(fpkm.ave$YLP0 != 0),]),
              YLP2.5=rownames(fpkm.ave[which(fpkm.ave$YLP2.5 != 0),]),
              YRP0=rownames(fpkm.ave[which(fpkm.ave$YRP0 != 0),]),
              YRP2.5=rownames(fpkm.ave[which(fpkm.ave$YRP2.5 != 0),]))

venn(PiSet, zcolor = "style")

# pdf("pics/Venn.pdf")
# venn(PiSet)
# dev.off()

geneList <- read.table("~/R/Eutrema/PS/Genes")
colnames(geneList) <- c("Family", "Name", "ID")


# graph these genes
# btw these are my Pi/S genes lmao
# maybe I should put some lipid biosynthesis genes in here 
library(ggplot2)
library(reshape2)

gene.fpkm <- fpkm[which(rownames(fpkm) %in% geneList[,3]),]
gene.fpkm <- data.frame(ID=rownames(gene.fpkm), gene.fpkm)
gene.fpkm <- merge(geneList, gene.fpkm, by = "ID")
gene.fpkm <- gene.fpkm %>% arrange(Name)

gene.fp.melt <- melt(gene.fpkm, id.vars = c("ID", "Family",  "Name"))

{gene.plot <- ggplot(gene.fp.melt, aes(x = Name, y = log2(value + 1), fill = factor(variable))) +
                    geom_bar(stat = "identity", position = "dodge") +
                    theme(axis.text.x = element_text(angle = 45, hjust = 1))
                gene.plot}
# ggsave("pics/pigenes.pdf", dpi = 350, height = 6, width = 12)

# separate treatment and tissue
# Ptreat <- factor(rep(c("0", "2.5", "0", "2.5"), each = 2))
# tissue <- factor(rep(c("leaf", "root"), each = 4))
# sample <- c("YLP0-1",
#             "YLP0-2",
#             "YLP2.5-1",
#             "YLP2.5-2",
#             "YRP0-1",
#             "YRP0-2",
#             "YRP2.5-1",
#             "YRP2.5-2")
# 
# metadata <- data.frame(sampleNO = sample,
#                        Ptreat = Ptreat,
#                        tissue = tissue,
#                        libraryName = sample)
# 
# metadata <- mutate(metadata, countFile = paste0(metadata$libraryName, "/", metadata$libraryName,
#                                                 "QC.geneCounts.formatted.for.DESeq.txt"))
# 
# sampleTable <- data.frame(sampleName = metadata$libraryName,
#                           fileName = metadata$countFile,
#                           Ptreat = metadata$Ptreat,
#                           tissue = metadata$tissue,
#                           sampleNO = metadata$sampleNO)
# 
# 
# DESeq2Table <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
#                                           directory = "QoRTs",
#                                           design = ~ Ptreat + tissue)
# 
# keep <- rowSums(counts(DESeq2Table)) >= ncol(DESeq2Table)
# DESeq2Table <- DESeq2Table[keep,]
# 
# DESeq2Table$tissue <- relevel(DESeq2Table$tissue, ref = "leaf")
# DESeq2Table$Ptreat <- relevel(DESeq2Table$Ptreat, ref = "2.5")
# 
# dds <- DESeq(DESeq2Table)
# resultsNames(dds)


# Find the list of genes lmao

# Oh yeah, should also do a PCA

head(fpkm)
fpkm.log <- log2(fpkm + 1)
pca<- prcomp(fpkm.log, center=TRUE, scale=TRUE)  # PCA with centering and scaling
pca$rotation  # The loadings are here

sdev <- pca$sdev

screeplot(pca)
#log.expr <- log2(expr+1)

plot(pca, type = "l")
summary(pca)
exprVals<-data.frame(pca$x)
sampleVals<-data.frame(pca$rotation)
#rv <- rowVars(as.matrix(log.expr))
#names(rv) <- rownames(expr)
#select <- order(rv, decreasing = TRUE)[1:2000]
#pca <- prcomp(log.expr[select,], .scale=TRUE, center=TRUE)

dim(exprVals)
dim(sampleVals)

# graph the pca
library(ggrepel)
library(shadowtext)
library(extrafont)


samples <- c("YLP0-1", "YLP0-2", "YLP2.5-1", "YLP2.5-2", 
             "YRP0-1", "YRP0-2", "YRP2.5-1", "YRP2.5-2")
coords <- data.frame(X=rep(0, 8), Y=rep(0, 8),sampleVals, Samples = samples)
# coords$Treatment <- factor(coords$Samples, c("ps", "pS", "Ps", "PS"))

# extract the precent PC contribution
pcaSum <- as.data.frame(summary(pca)$importance)

set.seed(51)
# PC1 vs PC2
{fpkm.plot <- ggplot(exprVals, aes_string("PC1", "PC2")) +
    geom_point(shape=19, alpha=0.3) +
    #     geom_point(data = exprVals[which(toupper(rownames(exprVals)) %in% toupper(allGlnc)) ,], 
    #                colour="#9fcc2e", size=2, alpha = 0.7) +    
    geom_segment(data=coords, aes(x=X, y=Y, xend=PC1, yend=PC2, colour=Samples), 
                 arrow=arrow(length = unit(0.3, "cm"), angle = 45), size = 1.5) +
#   scale_colour_manual(values=c( "#febfcb", "gold", "#f91301", "#fd8a19" ), name = "") +
    xlab(paste0("PC1 ", pcaSum$PC1[2]*100,"%")) + ylab(paste0("PC2 ",pcaSum$PC2[2]*100,"%")) +
    coord_cartesian(xlim=c(-2,6), ylim=c(-1.5,1.5))+
theme(axis.title=element_text(size=15), legend.text=element_text(size=12),
axis.text=element_text(size=12))
fpkm.plot}

{fpkm.plot23 <- ggplot(exprVals, aes_string("PC2", "PC3")) +
    geom_point(shape=19, alpha=0.3) +
    #     geom_point(data = exprVals[which(toupper(rownames(exprVals)) %in% toupper(allGlnc)) ,], 
    #                colour="#9fcc2e", size=2, alpha = 0.7) +    
    geom_segment(data=coords, aes(x=X, y=Y, xend=PC2, yend=PC3, colour=Samples), 
                 arrow=arrow(length = unit(0.3, "cm"), angle = 45), size = 1.5) +
#   scale_colour_manual(values=c( "#febfcb", "gold", "#f91301", "#fd8a19" ), name = "") +
    xlab(paste0("PC2 ", pcaSum$PC2[2]*100,"%")) + ylab(paste0("PC3 ",pcaSum$PC3[2]*100,"%")) +
    coord_cartesian(xlim=c(-2,6), ylim=c(-1.5,1.5))+
theme(axis.title=element_text(size=15), legend.text=element_text(size=12),
axis.text=element_text(size=12))
fpkm.plot23}

## using fpkm.ave
head(fpkm.ave)
fpkm.log <- log2(fpkm.ave + 1)
pca<- prcomp(fpkm.log, center=TRUE, scale=TRUE)  # PCA with centering and scaling
pca$rotation  # The loadings are here

sdev <- pca$sdev

screeplot(pca)
#log.expr <- log2(expr+1)

plot(pca, type = "l")
summary(pca)
exprVals<-data.frame(pca$x)
sampleVals<-data.frame(pca$rotation)
#rv <- rowVars(as.matrix(log.expr))
#names(rv) <- rownames(expr)
#select <- order(rv, decreasing = TRUE)[1:2000]
#pca <- prcomp(log.expr[select,], .scale=TRUE, center=TRUE)

dim(exprVals)
dim(sampleVals)

# graph the pca
library(ggrepel)
library(shadowtext)
library(extrafont)


samples <- c("YLP0", "YLP2.5",
             "YRP0", "YRP2.5")
coords <- data.frame(X=rep(0, 4), Y=rep(0, 4),sampleVals, Samples = samples)
# coords$Treatment <- factor(coords$Samples, c("ps", "pS", "Ps", "PS"))

# extract the precent PC contribution
pcaSum <- as.data.frame(summary(pca)$importance)

set.seed(51)
# PC1 vs PC2
# {fpkm.plot <- ggplot(exprVals, aes_string("PC1", "PC2")) +
#     geom_point(shape=19, alpha=0.3) +
# geom_point(data = exprVals[which(toupper(rownames(exprVals)) %in% toupper(allGlnc)) ,], 
#            colour="#9fcc2e", size=2, alpha = 0.7) +    
#     geom_segment(data=coords, aes(x=X, y=Y, xend=PC1, yend=PC2, colour=Samples), 
#                  arrow=arrow(length = unit(0.3, "cm"), angle = 45), size = 1.5) +
# scale_colour_manual(values=c( "#febfcb", "gold", "#f91301", "#fd8a19" ), name = "") +
#     xlab(paste0("PC1 ", pcaSum$PC1[2]*100,"%")) + ylab(paste0("PC2 ",pcaSum$PC2[2]*100,"%")) +
#     coord_cartesian(xlim=c(-2,6), ylim=c(-1.5,1.5))+
# theme(axis.title=element_text(size=15), legend.text=element_text(size=12),
# axis.text=element_text(size=12))
# fpkm.plot}

# DEG names
R0R2.up
R0R2.dn
L2L0.up
L2L0.dn

DEGsnames <- rbind(data.frame(ID=rownames(R0R2.up), Treat=rep("Root Up", length(rownames(R0R2.up)))),
                   data.frame(ID=rownames(R0R2.dn), Treat=rep("Root Down", length(rownames(R0R2.dn)))),
                   data.frame(ID=rownames(L2L0.dn), Treat=rep("Leaf Up", length(rownames(L2L0.dn)))),
                   data.frame(ID=rownames(L2L0.dn), Treat=rep("Leaf Down", length(rownames(L2L0.dn)))))
# DEGsnames$Treat <- factor(DEGsnames$Treat, levels = unique(DEGsnames$Treat))

DEGsexprVals <- merge(DEGsnames, exprVals, by.x = "ID", by.y = 0)
# PC2 vs PC3: separates nicely
{fpkm.plot23 <- ggplot(exprVals, aes_string("PC2", "PC3")) +
    geom_point(shape=19, alpha=0.3) +
    #     geom_point(data = DEGsexprVals, aes("PC2", "PC3", color = DEGsexprVals$Treat)) +
    geom_point(data = exprVals[which(rownames(exprVals) %in% rownames(L2L0.up)) ,], 
               colour="#39db64", size=2, alpha = 0.7) +
geom_point(data = exprVals[which(rownames(exprVals) %in% rownames(L2L0.dn)) ,], 
           colour="#3579cc", size=2, alpha = 0.7) +
geom_point(data = exprVals[which(rownames(exprVals) %in% rownames(R0R2.up)) ,], 
           colour="#e8320e", size=2, alpha = 0.7) +
geom_point(data = exprVals[which(rownames(exprVals) %in% rownames(R0R2.dn)) ,], 
           colour="#e30b5a", size=2, alpha = 0.7) +
    #     scale_fill_manual(values = c("#39db64","#3579cc","#e8320e","#e30b5a"), name = "") +
    #     geom_point(data = exprVals[which(rownames(exprVals) %in% DEGsnames) ,], 
    #                colour="#9fcc2e", size=2, alpha = 0.7) +    
    geom_segment(data=coords, aes(x=X, y=Y, xend=PC2, yend=PC3, colour=Samples), 
                 arrow=arrow(length = unit(0.3, "cm"), angle = 45), size = 1.5) +
scale_colour_manual(values=c( "#febfcb", "gold", "#f91301", "#fd8a19" ), name = "") +
    xlab(paste0("PC2 ", pcaSum$PC2[2]*100,"%")) + ylab(paste0("PC3 ",pcaSum$PC3[2]*100,"%")) +
    coord_cartesian(xlim=c(-4,4), ylim=c(-1.5,1.5))+
theme(axis.title=element_text(size=15), legend.text=element_text(size=12),
axis.text=element_text(size=12))
fpkm.plot23}

ggsave("pics/fpkm-ave23.pdf", fpkm.plot23)

#######################################################################################
expr_annot <- read.csv("~/R/Eutrema/PS/FPKMS_LengthAdjusted.csv", row.names=1, header = T)
relevant_expr <- expr_annot[,c(1:5, 18,19,21,22)]

annot <- read.csv("~/Eutrema/FPKM/2018-10-15-drought_geneFPKM.csv", header = T)

lfcList <- list(Rup=R0R2.up, Rdn=R0R2.dn, Lup=L2L0.up, Ldn=L2L0.dn)

# remove duplicated 
annot <- annot %>% subset(!duplicated(gene_id))
rownames(annot) <- annot$gene_id
annot <- annot[,2:ncol(annot)]
annot <- annot[,c(1:8, 52:54, 56, 57)]

library(openxlsx)
library(xlsx)

# to handle GC overhead limit exceeded
# options(java.parameters = "- Xmx1024m")
 for (i in 1:length(lfcList)) {
     if ( i == 1) {
         app <- F
     } else {
         app <- T
     }
     lfcSet <- as.data.frame(lfcList[[i]])
     lfcSet$ID <- rownames(lfcSet)
     annotSubset <- annot[which(rownames(annot) %in% rownames(lfcSet)),]
     annotSubset$ID <- rownames(annotSubset)
     toWrite <- merge(lfcSet, annotSubset, by = "ID")
     write.xlsx2(toWrite, file = "DEGannotations.xlsx", 
                 append = app, sheetName = names(lfcList)[i], 
                 row.names = T, col.names = T)
 }

####################################################################3
# Now it's time to do a PCA!

fpkms <- gene.fpkm[,4:11]
PCs <- prcomp(log2(fpkms + 1), center = T, scale = T)
PCload <- data.frame(PCs$rotation)
pcSummary <- as.data.frame(summary(PCs)$importance) 
PCload$condition <- rep(c("L0", "L2.5", "R0", "R2.5"), each = 2)
PC12 <- PCload %>% ggplot(aes(x = PC1, y = PC2, color = condition)) + geom_point() + xlim(c(-1,1)) + ylim(c(-1,1))+ labs(title = "PC1 vs PC2" , y = paste0("PC2 (", pcSummary$PC2[2] * 100, "%)", x = paste0("PC1 (", pcSummary$PC1[2] * 100, "%)")))
PC33 <- PCload %>% ggplot(aes(x = PC3, y = PC3, color = condition)) + geom_point() + xlim(c(-1,1)) + ylim(c(-1,1)) + labs(title = "PC2 vs PC3" , y = paste0("PC3 (", pcSummary$PC3[2] * 100, "%)", x = paste0("PC2 (", pcSummary$PC2[2] * 100, "%)"))) 
library(cowplot)
ggsave("pics/PiPCAs.pdf", plot_grid(PC12, PC23, nrow = 1), dpi = 250)
