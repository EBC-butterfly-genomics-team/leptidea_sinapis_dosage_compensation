
setwd("G:/Masterarbete/DEseq/For_DC")
library(LSD)
library(RColorBrewer)
library(spatialEco)

instar_volcano <- read.csv("LFC_Instar.txt", header=TRUE, sep = "\t")
FBG_instar <- read.csv("gene_names-FBG_Instar.txt", header = TRUE, sep = "\t")
FBG_instar$gene_name <- gsub(FBG_instar$gene_name, pattern = '-', replacement = '')
MBG_instar <- read.csv("gene_names-MBG_Instar.txt", header = TRUE, sep = "\t")
MBG_instar$gene_name <- gsub(MBG_instar$gene_name, pattern = '-', replacement = '')

pupa_volcano <- read.csv("LFC_Pupa.txt", header=TRUE, sep = "\t")
FBG_pupa <- read.csv("gene_names-FBG_Pupa.txt", header = TRUE, sep = "\t")
FBG_pupa$gene_name <- gsub(FBG_pupa$gene_name, pattern = '-', replacement = '')
MBG_pupa <- read.csv("gene_names-MBG_Pupa.txt", header = TRUE, sep = "\t")
MBG_pupa$gene_name <- gsub(MBG_pupa$gene_name, pattern = '-', replacement = '')

adult_volcano <- read.csv("LFC_Adult.txt", header=TRUE, sep = "\t")
FBG_adult <- read.csv("gene_names-FBG_Adult.txt", header = TRUE, sep = "\t")
FBG_adult$gene_name <- gsub(FBG_adult$gene_name, pattern = '-', replacement = '')
MBG_adult <- read.csv("gene_names-MBG_Adult.txt", header = TRUE, sep = "\t")
MBG_adult$gene_name <- gsub(MBG_adult$gene_name, pattern = '-', replacement = '')

###########################################

#par(pty="s", mfrow = c(1,3))

colscale <- brewer.pal(9, "Greys")

heatscatter(instar_volcano$log2FoldChange, -log10(instar_volcano$padj), main = "",
            pch = 20, cexplot = 0.8, nrcol = 5, grid = 1000, colpal = colscale[4:8],
            xlim = c(-3,3), ylab = "-log10 p-adj", cex.lab = 1.5, xlab = "log2 fold change")
axis(3, labels = c("Instar V"), lwd = 0, cex.axis = 2, at = 0, tck = 0, line = 1)
points(MBG_instar$log2FoldChange, -log10(MBG_instar$padj), pch=20, cex = 1, col="blue")
points(FBG_instar$log2FoldChange, -log10(FBG_instar$padj), pch=20, cex = 1, col="red")
text(FBG_instar$log2FoldChange, -log10(FBG_instar$padj), labels = FBG_instar$gene_name,
     cex = 0.9, pos = 2)
text(MBG_instar$log2FoldChange, -log10(MBG_instar$padj), labels = MBG_instar$gene_name,
     cex = 0.9, pos = 4)
#instar_V_smad3 <- subset(MBG_instar, MBG_instar$gene_name=="smad3")
#text(instar_V_smad3$log2FoldChange, -log10(instar_V_smad3$padj), labels = instar_V_smad3$gene_name, cex = 1, pos = 4)

legend("topleft", legend = c("FBG", "MBG"), y.intersp = 0.5, cex = 1, pch = 19, col = c("red","blue"))
segments(-1.5, 0, -1.5, 110, lwd = 1, lty = 2)
segments(1.5, 0, 1.5, 110, lwd = 1, lty = 2)
segments(-20, -log10(0.05), 20, -log10(0.05), lwd = 1, lty = 2)

##########################################

heatscatter(pupa_volcano$log2FoldChange, -log10(pupa_volcano$padj), main = "", xlim = c(-10,10), pch = 20,
            cexplot = 0.8, nrcol = 5, grid = 1000, colpal = colscale[4:8],
            xlab = "log2 fold change", cex.lab = 1.5, ylab = "-log10 p-adj")
axis(3, labels = c("Pupa"), lwd = 0, cex.axis = 2, at = 0, tck = 0, line = 1)
points(MBG_pupa$log2FoldChange, -log10(MBG_pupa$padj), pch=20, cex = 1, col="blue")
points(FBG_pupa$log2FoldChange, -log10(FBG_pupa$padj), pch=20, cex = 1, col="red")
text(FBG_pupa$log2FoldChange, -log10(FBG_pupa$padj), labels = FBG_pupa$gene_name,
     cex = 0.9, pos = 2)
text(MBG_pupa$log2FoldChange, -log10(MBG_pupa$padj), labels = MBG_pupa$gene_name,
     cex = 0.9, pos = 4)
#pupa_smad3 <- subset(MBG_pupa, MBG_pupa$gene_name=="smad3")
#text(pupa_smad3$log2FoldChange, -log10(pupa_smad3$padj), labels = pupa_smad3$gene_name, cex = 1, pos = 4)

legend("topleft", legend = c("FBG", "MBG"), y.intersp = 0.5, cex = 1, pch = 19, col = c("red","blue"))
segments(-1.5, 0, -1.5, 110, lwd = 1, lty = 2)
segments(1.5, 0, 1.5, 110, lwd = 1, lty = 2)
segments(-20, -log10(0.05), 20, -log10(0.05), lwd = 1, lty = 2)

#########################################

heatscatter(adult_volcano$log2FoldChange, -log10(adult_volcano$padj), main = "", pch = 20,
            cexplot = 0.8, nrcol = 5, grid = 1000, colpal = colscale[4:8], xlim = c(-20,20),
            ylim = c(0,124), cex.lab = 1.5, xlab = "log2 fold change", ylab = "-log10 p-adj")
axis(3, labels = c("Adult"), lwd = 0, cex.axis = 2, at = 0, tck = 0, line = 1)
points(MBG_adult$log2FoldChange, -log10(MBG_adult$padj), pch=20, col="blue", cex = 1)
points(FBG_adult$log2FoldChange, -log10(FBG_adult$padj), pch=20, col="red", cex = 1)
text(FBG_adult$log2FoldChange, -log10(FBG_adult$padj), labels = FBG_adult$gene_name,
     cex = 0.9, pos = 2)
text(MBG_adult$log2FoldChange, -log10(MBG_adult$padj), labels = MBG_adult$gene_name,
     cex = 0.9, pos = 4)
#adult_smad3 <- subset(MBG_adult, MBG_adult$gene_name=="smad3")
#text(adult_smad3$log2FoldChange, -log10(adult_smad3$padj), labels = adult_smad3$gene_name, cex = 1, pos = 4)

legend("topleft", legend = c("FBG", "MBG"), y.intersp = 0.5, cex = 1, pch = 19, col = c("red","blue"))
segments(-1.5, 0, -1.5, 130, lwd = 1, lty = 2)
segments(1.5, 0, 1.5, 130, lwd = 1, lty = 2)
segments(-40, -log10(0.05), 40, -log10(0.05), lwd = 1, lty = 2)

#legend(0,-10, legend = "MBG", cex = 1.5, pch = 19, col = "blue", xpd = TRUE, bty = "n", x.intersp = 0.5)
#legend(-10,-10, legend = "FBG", cex = 1.5, pch = 19, col = "red", xpd = TRUE, bty = "n", x.intersp = 0.5)