library(DESeq2) #DEG analysis
library(dplyr) #data manipulation

wd <- "~/scratch/Pi_reads/"
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

# redos the nbinomWalTest with the new relevelled reference
ddsR0 <- dds
ddsR0$condition <- relevel(ddsR0$condition, ref = "R0")
ddsR0 <- nbinomWaldTest(ddsR0)
resultsNames(ddsR0)

R0_L0 <- results(ddsR0, contrast = c("condition", "L0", "R0"))
R0_L2.5 <- results(ddsR0, contrast = c("condition", "L2.5", "R0"))
R0_R2.5 <- results(ddsR0, contrast = c("condition", "R2.5", "R0"))

R0L0.up <- R0_L0[which(R0_L0[,"padj"] < 0.05 & R0_L0[,"log2FoldChange"] > 0),]
R0L0.dn <- R0_L0[which(R0_L0[,"padj"] < 0.05 & R0_L0[,"log2FoldChange"] < 0),]
R0L2.up <- R0_L2.5[which(R0_L2.5[,"padj"] < 0.05 & R0_L2.5[,"log2FoldChange"] > 0),]
R0L2.dn <- R0_L2.5[which(R0_L2.5[,"padj"] < 0.05 & R0_L2.5[,"log2FoldChange"] < 0),]
R0R2.up <- R0_R2.5[which(R0_R2.5[,"padj"] < 0.05 & R0_R2.5[,"log2FoldChange"] > 0),]
R0R2.dn <- R0_R2.5[which(R0_R2.5[,"padj"] < 0.05 & R0_R2.5[,"log2FoldChange"] < 0),]

nrow(R0L0.up)
nrow(R0L0.dn)
nrow(R0L2.up)
nrow(R0L2.dn)
nrow(R0R2.up)
nrow(R0R2.dn)

L0_R2.5 <- results(dds, contrast = c("condition", "R2.5", "L0"))

L0R2.up <- L0_R2.5[which(L0_R2.5[,"padj"] < 0.05 & L0_R2.5[,"log2FoldChange"] > 0),]
L0R2.dn <- L0_R2.5[which(L0_R2.5[,"padj"] < 0.05 & L0_R2.5[,"log2FoldChange"] < 0),]
nrow(L0R2.up)
nrow(L0R2.dn)


# Venn diagrams
library(GenomicFeatures)               #GFF functions makeTxdbFromGff
library(venn)

# get median length of libraries

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

pdf("pics/Venn.pdf")
venn(PiSet)
dev.off()

geneList <- read.table("~/R/Eutrema/PS/Genes")
colnames(geneList) <- c("Name", "ID")


# graph these genes

library(ggplot2)
library(reshape2)

gene.fpkm <- fpkm[which(rownames(fpkm) %in% geneList[,2]),]
gene.fpkm <- data.frame(ID=rownames(gene.fpkm), gene.fpkm)
gene.fpkm <- merge(geneList, gene.fpkm, by = "ID")
gene.fpkm <- gene.fpkm %>% arrange(Name)

gene.fp.melt <- melt(gene.fpkm, id.vars = c("ID", "Name"))

{gene.plot <- ggplot(gene.fp.melt, aes(x = Name, y = value, fill = factor(variable))) +
                    geom_bar(stat = "identity", position = "dodge") +
                    theme(axis.text.x = element_text(angle = 45, hjust = 1))
                gene.plot}
ggsave("pics/pigenes.pdf", dpi = 350, height = 6, width = 12)

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



