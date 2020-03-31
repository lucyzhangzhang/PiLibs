#!/home/lucy/miniconda3/envs/renv/bin/Rscripts
library(dplyr)

wd <- "~/R/Eutrema/Pi"
setwd(wd)

# Arabidopsis dds
# Anames <- paste(rep(c("data/L2L0", "data/R2R0"), each = 2), c("up", "dn"), sep = ".")
Aname <- paste0("data/", rep(c(1,3), each = 2), "Day", rep(c("Root", "Shoot"), 2), "DEG.tab")
ADEGsList <- lapply(Aname, function(x) read.table(x))
names(ADEGsList) <- Aname
Atreat <- c("1DR", "1DS", "3DR", "3DS")

# Eutrema dds
library(xlsx)
library(readxl)
library(openxlsx)
Edata <- "~/R/Eutrema/PS/DEGsList.xlsx"
Enames <- read_excel(Edata)
Esheets <- excel_sheets(Edata)
EDegsList <- lapply(Esheets, read_excel, path = Edata)

names(EDegsList) <- Esheets

EDegsList <- lapply(EDegsList, function(x) {
                        x$Ara <- gsub("\\.[0-9]*", "", x$Ara)
                        return(x)
})

# DEGsCombined <- lapply(EDegsList, function(x) {
#         x$Overlap <- 
# })

test <- EDegsList$A_S_vs_s$Ara
test <- test[!is.na(test)]

inOrNot <- function(term) {
    stats <- c()
    if (!is.na(term)) {
    	for (i in 1:4) {
        	sig <- ADEGsList[[i]]
        	sig <- sig[which(sig$padj <= 0.05),]
        	sigU <- sig[which(sig$log2FoldChange > 0),]
        	sigD <- sig[which(sig$log2FoldChange < 0),]
        	stat <- ifelse(term %in% rownames(sigU), paste0(Atreat[i], "Up"), 
        	               ifelse(term %in% rownames(sigD), paste0(Atreat[i], "Down"), list(NULL)))
                stats <- c(stats, unlist(stat))
    	}
    } else {
    stats <- NULL
    }
    res <- paste(stats, collapse = ":")
    return(res)

}

annotation <- read.csv("~/Eutrema/FPKM/2018-10-15-drought_geneFPKM.csv", header = T)
annotation <- annotation[!duplicated(annotation$gene_id),]
PSlib <- read.table("~/R/Eutrema/PS/FPKMS_LengthAdjusted2.tab", header = T)
write.csv(annotation, file = "~/Eutrema/FPKM/2018-drough-annot2.csv", row.names = F)
getAnnot <- function(term) {
    getA <- annotation[which(annotation$gene_id == term), 58]
    return(getA)
}

res <- lapply(EDegsList, function(x) {
           x$annot <- sapply(x$Eutr, getAnnot)
           x$overlap <- sapply(x$Ara, inOrNot)
           return(x[,2:6])
                       })

# write the data into excel for Elizabeth
# Can't do apply because of the append option
for (i in 1:length(res)) {
    if (i == 1) {
        app <- F
    } else {
        app <- T
    }
    # need to change it to a data.frame because tibbles doesn't support row.names
    # row.names/col.names = F option breaks with tibbles
    write.xlsx2(data.frame(res[[i]]), file = "data/AthvsEutr.xlsx", 
                sheetName = names(res)[i], append = app, row.names = F)
} 

m395 <- annotation[which(annotation$gene_id == "DROUGHT.25549"),][c(10:51)]
library(reshape2)
library(ggplot2)
library(dplyr)
mm395 <- melt(m395)
# plot(mm395)
{m395plot <- mm395 %>% ggplot(., aes(x = variable, y = value, fill = variable)) + geom_bar(stat = "identity", position = "dodge") + labs(title = "Expression of mir395e/DROUGHT.25549 in drought lib", x = "Library name", y = "FPKM") + theme(axis.text.x = element_text(angle = -45, hjust = -0.01), legend.position = "none")}

ggsave("pics/testplot.pdf", m395plot, height = 7, width = 10)
# None lol
ps395 <- PSlib[which(PSlib$gene_id == "DROUGHT.25549"),]

# Another mir824/Thhalv10027603m.g
m824 <- annotation[which(annotation$gene_id == "Thhalv10027603m.g"),][c(10:51)]
mm824 <- melt(m824)
# plot(mm824)
{m824plot <- mm824 %>% ggplot(., aes(x = variable, y = value, fill = variable)) + geom_bar(stat = "identity", position = "dodge") + labs(title = "Expression of mir824/Thhalv10027603m.g in drought lib", x = "Library name", y = "FPKM") + theme(axis.text.x = element_text(angle = -45, hjust = -0.01), legend.position = "none")} # 

ggsave("pics/m824.pdf", m824plot, height = 7, width = 10)

ps824 <- PSlib[which(PSlib$Gene == "Thhalv10027603m.g"),]
psm824 <- melt(ps824[2:length(ps824)])
psm824$Treat <- as.factor(rep(c("ps", "Ps", "pS", "PS"), each = 3))
{ps824plot <- psm824 %>% ggplot(., aes(x = variable, y = value, fill = Treat)) + 
    labs(title = "Expression of mir824/Thhalv10027603m.g in PS lib", x = "Library name", y = "FPKM") + 
    geom_bar(position = "dodge", stat = "identity") +
    scale_fill_manual(values=c( "#febfcb", "gold", "#f91301", "#fd8a19" ))
}
ps824plot
library(grid)
library(gridExtra)
library(cowplot)
combinedplot <- plot_grid(m824plot, m395plot, ps824plot, ncol = 2, alight = "v")
combinedplot
ggsave("pics/combined.pdf", combinedplot, height = 14, width = 18)

# Comparing the Arabidopsis DEGs to the Pi libraries
# They are all padj significant
Pname <- paste0("data/",rep(c("L2L0", "R2R0"), each = 2), c(".up", ".dn"))
pPname <- paste0(rep(c("L2L0", "R2R0"), each = 2), c(".up", ".dn"))
PDEGsList <- lapply(Pname, function(x) read.table(x, row.names = 1))
names(PDEGsList) <- pPname

araConvert <- read.table("~/Eutrema/FPKM/eutremaToArabidopsis.names")
colnames(araConvert) <- c("Eutr", "Ara")
araConvert$Ara <- gsub("\\.[0-9]*", "", araConvert$Ara)

PDEGConvert <- lapply(PDEGsList, function(x) {
                          x$Expr <- ifelse(x[,2] > 0, "Up", "Down")
                          x <- merge(x, araConvert, by.x = 0, by.y = "Eutr", all.x = T)
                          x$Ara <- gsub("\\.[0-9]*", "", x$Ara)
                          x$annot <- sapply(x$Row.names, getAnnot)
                          x$overlap <- sapply(x$Ara, inOrNot)
                          rownames(x) <- x$Row.names
                          return(x[,c(9,8, 10,11)])
                })

for (i in 1:length(PDEGConvert)) {
    if (i == 1) {
        app <- F
    } else {
        app <- T
    }
    write.xlsx2(PDEGConvert[[i]], file = "data/PiLibComparison.xlsx",  sheetName = pPname[i], append = app, rowNames = T, quote = F)

}

# read the data coming from the tables


morcuende.1 <- read.xlsx2("data/morcuende2006.xlsx",  1, header = F)
morcuende.2 <- read.xlsx2("data/morcuende2006.xlsx",  2, header = F)
morcuende.raw <- rbind(morcuende.1, morcuende.2)
colnames(morcuende.raw) <- c("Affymetrix", "Ara", "Description", "Other nutrients")

library(stringr)
morcuende.trimwhite <- as.data.frame(apply(morcuende.raw, 2, str_trim, side = "both"))
morcuende.trimwhite$Ara <- toupper(morcuende.trimwhite$Ara)

misson.raw <- read.xlsx2("data/misson2005.xlsx", 1, header = F)
misson.raw <- misson.raw[,1:4]
colnames(misson.raw) <- c("Eutr", "Ara", "Name", "Description")
misson.raw$Name <- gsub(".*;", "", misson.raw$Name)
misson.raw$Ara <- gsub("\\.[0-9]*", "", misson.raw$Ara)

# pnas.raw <- read.csv("data/pnas_0505266102_05266Table5.csv", header = T)

woo2012 <- read.xlsx2("data/woo2012.xlsx", 6, header = T)

PSRs <- data.frame(Ara=c(as.character(woo2012$Ara), morcuende.trimwhite$Ara, misson.raw$Ara), stringsAsFactors = F)
PSRs <- unique(PSRs)
PSRs <- merge(araConvert, PSRs, by = "Ara")
PSRs$annot <- sapply(PSRs$Eutr, getAnnot)

PHTs <- read.csv("~/Eutrema/FPKM/PHTs.csv", header = F)
PHTs <- PHTs[,c(2,1,4)]
colnames(PHTs) <- c("Ara", "Eutr", "annot")
PHTs$Ara <- gsub("\\.[0-9]*", "", PHTs$Ara)

PSRs <- rbind(PHTs, PSRs)

PSR.big <- lapply(PDEGsList, function(x) {
                          x$Expr <- ifelse(x[,2] > 0, "Up", "Down")
                          x <- merge(x, PSRs, by.x = 0, by.y = "Eutr")
                          #                           x <- merge(x, PSRs, by.x = 0, by.y = "Eutr", all.x = T)
                          rownames(x) <- x$Row.names
                          return(x[,c(9, 8, 10)])
                })
lapply(PSR.big, tail)

for (i in 1:length(PSR.big)) {
    if (i == 1) {
        app <- F
    } else {
        app <- T
    }
    write.xlsx2(PSR.big[[i]], file = "data/KeyPSTRs.xlsx",  sheetName = pPname[i], append = app, rowNames = T, quote = F)

}

# same as the first one except remove Thhavls that are not also expressed in Ara
PDEGConvert2 <- lapply(PDEGsList, function(x) {
                          x$Expr <- ifelse(x[,2] > 0, "Up", "Down")
                          x <- merge(x, araConvert, by.x = 0, by.y = "Eutr")
                          #                           x <- merge(x, araConvert, by.x = 0, by.y = "Eutr", all.x = T)
                          #                           x$Ara <- gsub("\\.[0-9]*", "", x$Ara)
                          x$annot <- sapply(x$Row.names, getAnnot)
                          x$overlap <- sapply(x$Ara, inOrNot)
                          rownames(x) <- x$Row.names
                          return(x[,c(9,8, 10,11)])
                })

PDEGConvert2 <- lapply(PDEGConvert2, function(x) {
                           x <- x[which(x$overlap != ""),]
                           return(x)
                })

for (i in 1:length(PDEGConvert2)) {
    if (i == 1) {
        app <- F
    } else {
        app <- T
    }
    write.xlsx2(PDEGConvert2[[i]], file = "data/PiLibComparison2.xlsx",  sheetName = pPname[i], append = app, rowNames = T, quote = F)

}
