# Differential gene expression using DEseq2

#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
library("DESeq2")
library("IHW")
setwd("G:/Masterarbete/DEseq")

#Load gene count matrix and labels ####
countData <- as.matrix(read.csv("gene_count_matrix.csv", row.names="gene_id"))
colData <- read.csv("sample_info.txt", sep="\t", row.names=1)
#Note: The PHENO_DATA file contains information on each sample, e.g., sex or population.
#The exact way to import this depends on the format of the file.

#From deseq2 manual:
#...for some datasets, exploratory data analysis (EDA) plots could reveal that one
#or more groups has much higher within-group variability than the others. A simulated
#example of such a set of samples is shown below. This is case where, by comparing 
#groups A and B separately - subsetting a DESeqDataSet to only samples from those two 
#groups and then running DESeq on this subset - will be more sensitive than a model 
#including all samples together. It should be noted that such an extreme range of 
#within-group variability is not common, although it could arise if certain treatments
#produce an extreme reaction (e.g. cell death). Again, this can be easily detected
#from the EDA plots such as PCA described in this vignette.

#Make separate data sets for each developmental stage ####
Instar_countdata <- subset(countData, select = c("P5052_202_S48", "P5052_210_S51", "P5052_218_S55", 
                                                 "P5052_203_S49", "P5052_211_S52", "P5052_219_S56"))
Instar_coldata <- as.data.frame(colData[c(1, 5, 9, 2, 6, 10), ])

Pupa_countdata <- subset(countData, select = c("P5052_204_S34", "P5052_212_S36", "P5052_220_S57",
                                               "P5052_205_S35", "P5052_213_S53", "P5052_221_S58"))
Pupa_coldata <- as.data.frame(colData[c(3, 7, 11, 4, 8, 12)  , ])

Adult_countdata <- subset(countData, select = c("P5052_233_S60", "P5052_234_S61", "P5052_235_S62",
                                                "P5052_226_S59", "P5052_227_S38", "P7553_325_S10"))
Adult_coldata <- as.data.frame(colData[c(15, 16, 17, 13, 14, 18), ])

Instar_coldata <- Instar_coldata[ , 2:3]
Pupa_coldata <- Pupa_coldata[ , 2:3]
Adult_coldata <- Adult_coldata[ , 2:3]

#Pre-filtering I: Only keep genes with >1 samples with non-zero read count ####
Instar_countdata <- Instar_countdata[rowSums(Instar_countdata > 0) >= 2 , ]
Pupa_countdata <- Pupa_countdata[rowSums(Pupa_countdata > 0) >= 2 , ]
Adult_countdata <- Adult_countdata[rowSums(Adult_countdata > 0) >= 2 , ]

#Check all sample IDs in colData are also in CountData and match their orders ####
all(rownames(colData) %in% colnames(countData))
countData <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countData))

all(rownames(Instar_coldata) %in% colnames(Instar_countdata))
Instar_countdata <- Instar_countdata[, rownames(Instar_coldata)]
all(rownames(Instar_coldata) == colnames(Instar_countdata))

all(rownames(Pupa_coldata) %in% colnames(Pupa_countdata))
Pupa_countdata <- Pupa_countdata[, rownames(Pupa_coldata)]
all(rownames(Pupa_coldata) == colnames(Pupa_countdata))

all(rownames(Adult_coldata) %in% colnames(Adult_countdata))
Adult_countdata <- Adult_countdata[, rownames(Adult_coldata)]
all(rownames(Adult_coldata) == colnames(Adult_countdata))


#Create a DESeqDataSet from count matrix and labels ####
ddsInstar <- DESeqDataSetFromMatrix(countData = Instar_countdata,
                                    colData = Instar_coldata, design = ~Family + Sex)
ddsPupa <- DESeqDataSetFromMatrix(countData = Pupa_countdata,
                                  colData = Pupa_coldata, design = ~Family + Sex)
#Adult have no family confounding effects, see pca plot. Don't adjust for family...
ddsAdult <- DESeqDataSetFromMatrix(countData = Adult_countdata,
                                   colData = Adult_coldata, design = ~Sex)

#Pre-filtering II: Only keep rows with baseMean >1
ddsInstar <- ddsInstar[rowMeans(counts(ddsInstar)) > 1, ]
ddsPupa <- ddsPupa[rowMeans(counts(ddsPupa)) > 1, ]
ddsAdult <- ddsAdult[rowMeans(counts(ddsAdult)) > 1, ]


# PCA ####
library(ggplot2)
library(cowplot)

vsdInstar <- vst(ddsInstar)

pcaDataInstar <- plotPCA(vsdInstar, intgroup=c("Family", "Sex"), returnData=TRUE)
percentVarInstar <- round(100 * attr(pcaDataInstar, "percentVar"))
Instar_V <- ggplot(pcaDataInstar, aes(PC1, PC2, color=Sex, shape=Family)) +
  geom_point(size=2, stroke=1) +
  scale_color_manual(values = c("red", "blue")) +
  scale_shape_manual(values = c(1, 3, 8)) +
  xlab(paste0("PC1: ",percentVarInstar[1],"%")) +
  ylab(paste0("PC2: ",percentVarInstar[2],"%")) + 
  ggtitle("Instar V") +
  xlim(-120, 120) +
  ylim(-70, 70) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "plain"),
        panel.background = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(fill = NA), axis.text = element_text(size = 9),
        axis.title=element_text(size=10), legend.text=element_text(size=8),
        legend.title=element_text(size=10)) +
  coord_fixed()


vsdPupa <- vst(ddsPupa)

pcaDataPupa <- plotPCA(vsdPupa, intgroup=c("Family", "Sex"), returnData=TRUE)
percentVarPupa <- round(100 * attr(pcaDataPupa, "percentVar"))
Pupa <- ggplot(pcaDataPupa, aes(PC1, PC2, color=Sex, shape=Family)) +
  geom_point(size=2, stroke=1) +
  scale_color_manual(values = c("red", "blue")) +
  scale_shape_manual(values = c(1, 3, 8)) +
  xlab(paste0("PC1: ",percentVarPupa[1],"%")) +
  ylab(paste0("PC2: ",percentVarPupa[2],"%")) + 
  ggtitle("Pupa") +
  xlim(-120, 120) +
  ylim(-70, 70) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "plain"),
        panel.background = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(fill = NA), axis.text = element_text(size = 9),
        axis.title=element_text(size=10), legend.text=element_text(size=8),
        legend.title=element_text(size=10)) +
  coord_fixed()

vsdAdult <- vst(ddsAdult)

pcaDataAdult <- plotPCA(vsdAdult, intgroup=c("Family", "Sex"), returnData=TRUE)
percentVarAdult <- round(100 * attr(pcaDataAdult, "percentVar"))
Adult <- ggplot(pcaDataAdult, aes(PC1, PC2, colour=Sex, shape=Family)) +
  geom_point(size=2, stroke=1) +
  scale_color_manual(values = c("red", "blue")) +
  scale_shape_manual(values = c(1, 8, 0)) +
  xlab(paste0("PC1: ",percentVarAdult[1],"%")) +
  ylab(paste0("PC2: ",percentVarAdult[2],"%")) + 
  ggtitle("Adult") +
  xlim(-120, 120) +
  ylim(-70, 70) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "plain"),
        panel.background = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(fill = NA), axis.text = element_text(size = 9),
        axis.title=element_text(size=10), legend.text=element_text(size=8),
        legend.title=element_text(size=10)) +
  coord_fixed()

plot_grid(Instar_V, Pupa, Adult, ncol = 3)

#Run the default analysis for DESeq2 ####

#dds <- DESeq(dds)
ddsInstar <- DESeq(ddsInstar)
ddsPupa <- DESeq(ddsPupa)
ddsAdult <- DESeq(ddsAdult)

#Results ####

#Independent hypothesis weighting
resInstarIHW <- results(ddsInstar, filterFun = ihw, alpha = 0.05)
resPupaIHW <- results(ddsPupa, filterFun = ihw, alpha = 0.05)
resAdultIHW <- results(ddsAdult, filterFun = ihw, alpha = 0.05)

#Sort by adjusted p-value and display ####
(resInstarOrdered <- resInstarIHW[order(resInstarIHW$padj), ])
(resPupaOrdered <- resPupaIHW[order(resPupaIHW$padj), ])
(resAdultOrdered <- resAdultIHW[order(resAdultIHW$padj), ])

#summary(resInstar)
summary(resInstarIHW)
#summary(resPupa)
summary(resPupaIHW)
#summary(resAdult)
summary(resAdultIHW)

# Filter p-adj (FDR) < 0.05 ####
Filtered_P0.05_Instar <- subset(resInstarIHW, padj<0.05)
Filtered_P0.05_Pupa <- subset(resPupaIHW, padj<0.05)
Filtered_P0.05_Adult <- subset(resAdultIHW, padj<0.05)

summary(Filtered_P0.05_Instar)
summary(Filtered_P0.05_Pupa)
summary(Filtered_P0.05_Adult)

# Filter logFC > |1.5|) ####
Filtered_P0.05_LFC_1.5_Instar <- subset(Filtered_P0.05_Instar, log2FoldChange >1.5 | log2FoldChange < -1.5)
Filtered_P0.05_LFC_1.5_Pupa <- subset(Filtered_P0.05_Pupa, log2FoldChange >1.5 | log2FoldChange < -1.5)
Filtered_P0.05_LFC_1.5_Adult <- subset(Filtered_P0.05_Adult, log2FoldChange >1.5 | log2FoldChange < -1.5)

summary(Filtered_P0.05_LFC_1.5_Instar)
summary(Filtered_P0.05_LFC_1.5_Pupa)
summary(Filtered_P0.05_LFC_1.5_Adult)

#Filter baseMean > 10 ####

Filtered_P0.05_LFC_1.5_base_Mean_10_Instar <- subset(Filtered_P0.05_LFC_1.5_Instar, baseMean>10)
Filtered_P0.05_LFC_1.5_base_Mean_10_Pupa <- subset(Filtered_P0.05_LFC_1.5_Pupa, baseMean>10)
Filtered_P0.05_LFC_1.5_base_Mean_10_Adult <- subset(Filtered_P0.05_LFC_1.5_Adult, baseMean>10)

summary(Filtered_P0.05_LFC_1.5_base_Mean_10_Instar)
summary(Filtered_P0.05_LFC_1.5_base_Mean_10_Pupa)
summary(Filtered_P0.05_LFC_1.5_base_Mean_10_Adult)


# Select Male/Female biased genes ####
Filtered_P0.05_MBG_Instar <- subset(Filtered_P0.05_LFC_1.5_base_Mean_10_Instar, log2FoldChange >1.5)
Filtered_P0.05_MBG_Pupa <- subset(Filtered_P0.05_LFC_1.5_base_Mean_10_Pupa , log2FoldChange > 1.5)
Filtered_P0.05_MBG_Adult <- subset(Filtered_P0.05_LFC_1.5_base_Mean_10_Adult, log2FoldChange > 1.5)
summary(Filtered_P0.05_MBG_Instar)
summary(Filtered_P0.05_MBG_Pupa)
summary(Filtered_P0.05_MBG_Adult)

Filtered_P0.05_FBG_Instar <- subset(Filtered_P0.05_LFC_1.5_base_Mean_10_Instar, log2FoldChange < -1.5)
Filtered_P0.05_FBG_Pupa <- subset(Filtered_P0.05_LFC_1.5_base_Mean_10_Pupa, log2FoldChange < -1.5)
Filtered_P0.05_FBG_Adult <- subset(Filtered_P0.05_LFC_1.5_base_Mean_10_Adult, log2FoldChange < -1.5)
summary(Filtered_P0.05_FBG_Instar)
summary(Filtered_P0.05_FBG_Pupa)
summary(Filtered_P0.05_FBG_Adult)

# Filter out genes without expression
expressed_instar <- subset(resInstarIHW, baseMean >0)
expressed_pupa <- subset(resPupaIHW, baseMean >0)
expressed_adult <- subset(resAdultIHW, baseMean >0)

# Write results to file ####
write.table(Filtered_P0.05_LFC_1.5_base_Mean_10_Instar,
          file = "Filtered_LFC_1.5_baseMean_10_Instar.txt", sep = "\t", quote = FALSE)
write.table(Filtered_P0.05_LFC_1.5_base_Mean_10_Pupa, 
          file = "Filtered_LFC_1.5_baseMean_10_Pupa.txt", sep = "\t", quote = FALSE)
write.table(Filtered_P0.05_LFC_1.5_base_Mean_10_Adult, 
          file = "Filtered_LFC_1.5_baseMean_10_Adult.txt", sep = "\t", quote = FALSE)

write.table(Filtered_P0.05_MBG_Instar, 
          file = "MBG_Instar.txt", sep = "\t", quote = FALSE)
write.table(Filtered_P0.05_MBG_Pupa, 
          file = "MBG_Pupa.txt", sep = "\t", quote = FALSE)
write.table(Filtered_P0.05_MBG_Adult, 
          file = "MBG_Adult.txt", sep = "\t", quote = FALSE)

write.table(Filtered_P0.05_FBG_Instar, 
          file = "FBG_Instar.txt", sep = "\t", quote = FALSE)
write.table(Filtered_P0.05_FBG_Pupa, 
          file = "FBG_Pupa.txt", sep = "\t", quote = FALSE)
write.table(Filtered_P0.05_FBG_Adult, 
          file = "FBG_Adult.txt", sep = "\t", quote = FALSE)

#All LFC
write.table(expressed_instar, file = "LFC_Instar.txt", sep = "\t", quote = FALSE)
write.table(expressed_pupa, file = "LFC_Pupa.txt", sep = "\t", quote = FALSE)
write.table(expressed_adult, file = "LFC_Adult.txt", sep = "\t", quote = FALSE)
