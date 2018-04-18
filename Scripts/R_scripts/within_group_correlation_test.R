#Correlation test of gene expression as FPKM within each group

library("psych")

setwd("G:/Masterarbete/samples")

#######################################################################################

##### Instar V female ######

S49 <- read.delim("instar_V_female/S49.gene_abund.tab", sep = "\t", header = TRUE)
S52 <- read.delim("instar_V_female/S52.gene_abund.tab", sep = "\t", header = TRUE)
S56 <- read.delim("instar_V_female/S56.gene_abund.tab", sep = "\t", header = TRUE)

log_s49 <-log2(S49$FPKM)
log_s52 <-log2(S52$FPKM)
log_s56 <-log2(S56$FPKM)

rho_S49_S52 <- cor.test(log_s49, log_s52, method = "spearman", exact = FALSE)
rho_S49_S52 <- round(rho_S49_S52$estimate, digits = 2)
rho_S49_S56 <- cor.test(log_s49, log_s56, method = "spearman", exact = FALSE)
rho_S49_S56 <- round(rho_S49_S56$estimate, digits = 2)
rho_S52_S56 <- cor.test(log_s52, log_s56, method = "spearman", exact = FALSE)
rho_S52_S56 <- round(rho_S52_S56$estimate, digits = 2)

par(pty="s", mfrow = c(2,3))

plot(log_s49, log_s52, asp = TRUE, pch = 16, cex = 0.3,
     ylim = c(-8,18), xlim = c(-8,18), 
     title(main = "A", adj = 0, cex.main = 2), xlab = "", ylab = "")
axis(1, labels = "s49 log2 FPKM", cex.axis = 1.5, at = 5, line = 1.5, tck = 0) #x
axis(2, labels = "s52 log2 FPKM", cex.axis = 1.5, at = 5, line = 1.5, tck = 0) #y

text(13.8, -5.1, labels="0.90", cex = 1.2)
text(10,-5, "\\*r", vfont = c("sans serif", "bold italic"), cex=2)
text(11.5,-5, "=", vfont = c("sans serif", "plain"), cex = 2)

plot(log_s49, log_s56, asp = TRUE, pch = 16, cex = 0.3,
     ylim = c(-8,18), xlim = c(-8,18),
     title(main = "", adj = 0, cex.main = 2), xlab = "", ylab = "")
axis(1, labels = "s49 log2 FPKM", cex.axis = 1.5, at = 5, line = 1.5, tck = 0) #x
axis(2, labels = "s56 log2 FPKM", cex.axis = 1.5, at = 5, line = 1.5, tck = 0) #y

text(13.8, -5.1, labels=rho_S49_S56, cex = 1.2)
text(10,-5, "\\*r", vfont = c("sans serif", "bold italic"), cex=2)
text(11.5,-5, "=", vfont = c("sans serif", "plain"), cex = 2)

plot(log_s52, log_s56, asp = TRUE, pch = 16, cex = 0.3,
     ylim = c(-8,18), xlim = c(-8,18),
     title(main = "", adj = 0, cex.main = 2), xlab = "", ylab = "")
axis(1, labels = "s52 log2 FPKM", cex.axis = 1.5, at = 5, line = 1.5, tck = 0) #x
axis(2, labels = "s56 log2 FPKM", cex.axis = 1.5, at = 5, line = 1.5, tck = 0) #y

text(13.8, -5.1, labels="0.90", cex = 1.2)
text(10,-5, "\\*r", vfont = c("sans serif", "bold italic"), cex=2)
text(11.5,-5, "=", vfont = c("sans serif", "plain"), cex = 2)

##### Instar V Male ######

S48 <- read.delim("instar_V_male/S48.gene_abund.tab", sep = "\t", header = TRUE)
S51 <- read.delim("instar_V_male/S51.gene_abund.tab", sep = "\t", header = TRUE)
S55 <- read.delim("instar_V_male/S55.gene_abund.tab", sep = "\t", header = TRUE)

log_s48 <-log2(S48$FPKM)
log_s51 <-log2(S51$FPKM)
log_s55 <-log2(S55$FPKM)

rho_S48_S51 <- cor.test(log_s48, log_s51, method = "spearman", exact = FALSE)
rho_S48_S51 <- round(rho_S48_S51$estimate, digits = 2)
rho_S48_S55 <- cor.test(log_s48, log_s55, method = "spearman", exact = FALSE)
rho_S48_S55 <- round(rho_S48_S55$estimate, digits = 2)
rho_S51_S55 <- cor.test(log_s51, log_s55, method = "spearman", exact = FALSE)
rho_S51_S55 <- round(rho_S51_S55$estimate, digits = 2)

plot(log_s48, log_s51, asp = TRUE, pch = 16, cex = 0.3,
     ylim = c(-8,18), xlim = c(-8,18), 
     title(main = "", adj = 0, cex.main = 2), xlab = "", ylab = "")
axis(1, labels = "s48 log2 FPKM", cex.axis = 1.5, at = 5, line = 1.5, tck = 0) #x
axis(2, labels = "s51 log2 FPKM", cex.axis = 1.5, at = 5, line = 1.5, tck = 0) #y

text(13.8, -5.1, labels=rho_S48_S51, cex = 1.2)
text(10,-5, "\\*r", vfont = c("sans serif", "bold italic"), cex=2)
text(11.5,-5, "=", vfont = c("sans serif", "plain"), cex = 2)

plot(log_s48, log_s55, asp = TRUE, pch = 16, cex = 0.3,
     ylim = c(-8,18), xlim = c(-8,18),
     title(main = "", adj = 0, cex.main = 2), xlab = "", ylab = "")
axis(1, labels = "s48 log2 FPKM", cex.axis = 1.5, at = 5, line = 1.5, tck = 0) #x
axis(2, labels = "s55 log2 FPKM", cex.axis = 1.5, at = 5, line = 1.5, tck = 0) #y

text(13.8, -5.1, labels=rho_S48_S55, cex = 1.2)
text(10,-5, "\\*r", vfont = c("sans serif", "bold italic"), cex=2)
text(11.5,-5, "=", vfont = c("sans serif", "plain"), cex = 2)

plot(log_s51, log_s55, asp = TRUE, pch = 16, cex = 0.3,
     ylim = c(-8,18), xlim = c(-8,18),
     title(main = "", adj = 0, cex.main = 2), xlab = "", ylab = "")
axis(1, labels = "s51 log2 FPKM", cex.axis = 1.5, at = 5, line = 1.5, tck = 0) #x
axis(2, labels = "s55 log2 FPKM", cex.axis = 1.5, at = 5, line = 1.5, tck = 0) #y

text(13.8, -5.1, labels="0.90", cex = 1.2)
text(10,-5, "\\*r", vfont = c("sans serif", "bold italic"), cex=2)
text(11.5,-5, "=", vfont = c("sans serif", "plain"), cex = 2)

################################################################################

##### Pupa female ######

S35 <- read.delim("pupa_female/S35.gene_abund.tab", sep = "\t", header = TRUE)
S53 <- read.delim("pupa_female/S53.gene_abund.tab", sep = "\t", header = TRUE)
S58 <- read.delim("pupa_female/S58.gene_abund.tab", sep = "\t", header = TRUE)

log_s35 <-log2(S35$FPKM)
log_s53 <-log2(S53$FPKM)
log_s58 <-log2(S58$FPKM)

rho_S35_S53 <- cor.test(log_s35, log_s53, method = "spearman", exact = FALSE)
rho_S35_S53 <- round(rho_S35_S53$estimate, digits = 2)
rho_S35_S58 <- cor.test(log_s35, log_s58, method = "spearman", exact = FALSE)
rho_S35_S58 <- round(rho_S35_S58$estimate, digits = 2)
rho_S53_S58 <- cor.test(log_s53, log_s58, method = "spearman", exact = FALSE)
rho_S53_S58 <- round(rho_S53_S58$estimate, digits = 2)

par(pty="s", mfrow = c(2,3))

plot(log_s35, log_s53, asp = TRUE, pch = 16, cex = 0.3,
     ylim = c(-8,18), xlim = c(-8,18), 
     title(main = "B", adj = 0, cex.main = 2), xlab = "", ylab = "")
axis(1, labels = "s35 log2 FPKM", cex.axis = 1.5, at = 5, line = 1.5, tck = 0) #x
axis(2, labels = "s53 log2 FPKM", cex.axis = 1.5, at = 5, line = 1.5, tck = 0) #y

text(13.8, -5.1, labels=rho_S35_S53, cex = 1.2)
text(10,-5, "\\*r", vfont = c("sans serif", "bold italic"), cex=2)
text(11.5,-5, "=", vfont = c("sans serif", "plain"), cex = 2)

plot(log_s35, log_s58, asp = TRUE, pch = 16, cex = 0.3,
     ylim = c(-8,18), xlim = c(-8,18),
     title(main = "", adj = 0, cex.main = 2), xlab = "", ylab = "")
axis(1, labels = "s35 log2 FPKM", cex.axis = 1.5, at = 5, line = 1.5, tck = 0) #x
axis(2, labels = "s58 log2 FPKM", cex.axis = 1.5, at = 5, line = 1.5, tck = 0) #y

text(13.8, -5.1, labels=rho_S35_S58, cex = 1.2)
text(10,-5, "\\*r", vfont = c("sans serif", "bold italic"), cex=2)
text(11.5,-5, "=", vfont = c("sans serif", "plain"), cex = 2)

plot(log_s53, log_s58, asp = TRUE, pch = 16, cex = 0.3,
     ylim = c(-8,18), xlim = c(-8,18),
     title(main = "", adj = 0, cex.main = 2), xlab = "", ylab = "")
axis(1, labels = "s53 log2 FPKM", cex.axis = 1.5, at = 5, line = 1.5, tck = 0) #x
axis(2, labels = "s58 log2 FPKM", cex.axis = 1.5, at = 5, line = 1.5, tck = 0) #y

text(13.8, -5.1, labels=rho_S53_S58, cex = 1.2)
text(10,-5, "\\*r", vfont = c("sans serif", "bold italic"), cex=2)
text(11.5,-5, "=", vfont = c("sans serif", "plain"), cex = 2)

##### Pupa male ######

S34 <- read.delim("pupa_male/S34.gene_abund.tab", sep = "\t", header = TRUE)
S36 <- read.delim("pupa_male/S36.gene_abund.tab", sep = "\t", header = TRUE)
S57 <- read.delim("pupa_male/S57.gene_abund.tab", sep = "\t", header = TRUE)

log_s34 <-log2(S34$FPKM)
log_s36 <-log2(S36$FPKM)
log_s57 <-log2(S57$FPKM)

rho_S34_S36 <- cor.test(log_s34, log_s36, method = "spearman", exact = FALSE)
rho_S34_S36 <- round(rho_S34_S36$estimate, digits = 2)
rho_S34_S57 <- cor.test(log_s34, log_s57, method = "spearman", exact = FALSE)
rho_S34_S57 <- round(rho_S34_S57$estimate, digits = 2)
rho_S36_S57 <- cor.test(log_s36, log_s57, method = "spearman", exact = FALSE)
rho_S36_S57 <- round(rho_S36_S57$estimate, digits = 2)

plot(log_s34, log_s36, asp = TRUE, pch = 16, cex = 0.3,
     ylim = c(-8,18), xlim = c(-8,18), 
     title(main = "", adj = 0, cex.main = 2), xlab = "", ylab = "")
axis(1, labels = "s34 log2 FPKM", cex.axis = 1.5, at = 5, line = 1.5, tck = 0) #x
axis(2, labels = "s36 log2 FPKM", cex.axis = 1.5, at = 5, line = 1.5, tck = 0) #y

text(13.8, -5.1, labels=rho_S34_S36, cex = 1.2)
text(10,-5, "\\*r", vfont = c("sans serif", "bold italic"), cex=2)
text(11.5,-5, "=", vfont = c("sans serif", "plain"), cex = 2)

plot(log_s34, log_s57, asp = TRUE, pch = 16, cex = 0.3,
     ylim = c(-8,18), xlim = c(-8,18),
     title(main = "", adj = 0, cex.main = 2), xlab = "", ylab = "")
axis(1, labels = "s34 log2 FPKM", cex.axis = 1.5, at = 5, line = 1.5, tck = 0) #x
axis(2, labels = "s57 log2 FPKM", cex.axis = 1.5, at = 5, line = 1.5, tck = 0) #y

text(13.8, -5.1, labels=rho_S34_S57, cex = 1.2)
text(10,-5, "\\*r", vfont = c("sans serif", "bold italic"), cex=2)
text(11.5,-5, "=", vfont = c("sans serif", "plain"), cex = 2)

plot(log_s36, log_s57, asp = TRUE, pch = 16, cex = 0.3,
     ylim = c(-8,18), xlim = c(-8,18),
     title(main = "", adj = 0, cex.main = 2), xlab = "", ylab = "")
axis(1, labels = "s36 log2 FPKM", cex.axis = 1.5, at = 5, line = 1.5, tck = 0) #x
axis(2, labels = "s57 log2 FPKM", cex.axis = 1.5, at = 5, line = 1.5, tck = 0) #y

text(13.8, -5.1, labels=rho_S36_S57, cex = 1.2)
text(10,-5, "\\*r", vfont = c("sans serif", "bold italic"), cex=2)
text(11.5,-5, "=", vfont = c("sans serif", "plain"), cex = 2)

##### Adult female ######

S10 <- read.delim("adult_female/S10.gene_abund.tab", sep = "\t", header = TRUE)
S38 <- read.delim("adult_female/S38.gene_abund.tab", sep = "\t", header = TRUE)
S59 <- read.delim("adult_female/S59.gene_abund.tab", sep = "\t", header = TRUE)

log_s10 <-log2(S10$FPKM)
log_s38 <-log2(S38$FPKM)
log_s59 <-log2(S59$FPKM)

rho_S10_S38 <- cor.test(log_s10, log_s38, method = "spearman", exact = FALSE)
rho_S10_S38 <- round(rho_S10_S38$estimate, digits = 2)
rho_S10_S59 <- cor.test(log_s10, log_s59, method = "spearman", exact = FALSE)
rho_S10_S59 <- round(rho_S10_S59$estimate, digits = 2)
rho_S38_S59 <- cor.test(log_s38, log_s59, method = "spearman", exact = FALSE)
rho_S38_S59 <- round(rho_S38_S59$estimate, digits = 2)

par(pty="s", mfrow = c(2,3))

plot(log_s10, log_s38, asp = TRUE, pch = 16, cex = 0.3,
     ylim = c(-8,18), xlim = c(-8,18), 
     title(main = "C", adj = 0, cex.main = 2), xlab = "", ylab = "")
axis(1, labels = "s10 log2 FPKM", cex.axis = 1.5, at = 5, line = 1.5, tck = 0) #x
axis(2, labels = "s38 log2 FPKM", cex.axis = 1.5, at = 5, line = 1.5, tck = 0) #y

text(13.8, -5.1, labels=rho_S10_S38, cex = 1.2)
text(10,-5, "\\*r", vfont = c("sans serif", "bold italic"), cex=2)
text(11.5,-5, "=", vfont = c("sans serif", "plain"), cex = 2)

plot(log_s10, log_s59, asp = TRUE, pch = 16, cex = 0.3,
     ylim = c(-8,18), xlim = c(-8,18),
     title(main = "", adj = 0, cex.main = 2), xlab = "", ylab = "")
axis(1, labels = "s10 log2 FPKM", cex.axis = 1.5, at = 5, line = 1.5, tck = 0) #x
axis(2, labels = "s59 log2 FPKM", cex.axis = 1.5, at = 5, line = 1.5, tck = 0) #y

text(13.8, -5.1, labels=rho_S10_S59, cex = 1.2)
text(10,-5, "\\*r", vfont = c("sans serif", "bold italic"), cex=2)
text(11.5,-5, "=", vfont = c("sans serif", "plain"), cex = 2)

plot(log_s38, log_s59, asp = TRUE, pch = 16, cex = 0.3,
     ylim = c(-8,18), xlim = c(-8,18),
     title(main = "", adj = 0, cex.main = 2), xlab = "", ylab = "")
axis(1, labels = "s38 log2 FPKM", cex.axis = 1.5, at = 5, line = 1.5, tck = 0) #x
axis(2, labels = "s59 log2 FPKM", cex.axis = 1.5, at = 5, line = 1.5, tck = 0) #y

text(13.8, -5.1, labels=rho_S38_S59, cex = 1.2)
text(10,-5, "\\*r", vfont = c("sans serif", "bold italic"), cex=2)
text(11.5,-5, "=", vfont = c("sans serif", "plain"), cex = 2)

##### Adult male ######

S60 <- read.delim("adult_male/S60.gene_abund.tab", sep = "\t", header = TRUE)
S61 <- read.delim("adult_male/S61.gene_abund.tab", sep = "\t", header = TRUE)
S62 <- read.delim("adult_male/S62.gene_abund.tab", sep = "\t", header = TRUE)

log_s60 <-log2(S60$FPKM)
log_s61 <-log2(S61$FPKM)
log_s62 <-log2(S62$FPKM)

rho_S60_S61 <- cor.test(log_s60, log_s61, method = "spearman", exact = FALSE)
rho_S60_S61 <- round(rho_S60_S61$estimate, digits = 2)
rho_S60_S62 <- cor.test(log_s60, log_s62, method = "spearman", exact = FALSE)
rho_S60_S62 <- round(rho_S60_S62$estimate, digits = 2)
rho_S61_S62 <- cor.test(log_s61, log_s62, method = "spearman", exact = FALSE)
rho_S61_S62 <- round(rho_S61_S62$estimate, digits = 2)

plot(log_s60, log_s61, asp = TRUE, pch = 16, cex = 0.3,
     ylim = c(-8,18), xlim = c(-8,18), 
     title(main = "", adj = 0, cex.main = 2), xlab = "", ylab = "")
axis(1, labels = "s60 log2 FPKM", cex.axis = 1.5, at = 5, line = 1.5, tck = 0) #x
axis(2, labels = "s61 log2 FPKM", cex.axis = 1.5, at = 5, line = 1.5, tck = 0) #y

text(13.8, -5.1, labels=rho_S60_S61, cex = 1.2)
text(10,-5, "\\*r", vfont = c("sans serif", "bold italic"), cex=2)
text(11.5,-5, "=", vfont = c("sans serif", "plain"), cex = 2)

plot(log_s60, log_s62, asp = TRUE, pch = 16, cex = 0.3,
     ylim = c(-8,18), xlim = c(-8,18),
     title(main = "", adj = 0, cex.main = 2), xlab = "", ylab = "")
axis(1, labels = "s60 log2 FPKM", cex.axis = 1.5, at = 5, line = 1.5, tck = 0) #x
axis(2, labels = "s62 log2 FPKM", cex.axis = 1.5, at = 5, line = 1.5, tck = 0) #y

text(13.8, -5.1, labels=rho_S60_S62, cex = 1.2)
text(10,-5, "\\*r", vfont = c("sans serif", "bold italic"), cex=2)
text(11.5,-5, "=", vfont = c("sans serif", "plain"), cex = 2)

plot(log_s61, log_s62, asp = TRUE, pch = 16, cex = 0.3,
     ylim = c(-8,18), xlim = c(-8,18),
     title(main = "", adj = 0, cex.main = 2), xlab = "", ylab = "")
axis(1, labels = "s61 log2 FPKM", cex.axis = 1.5, at = 5, line = 1.5, tck = 0) #x
axis(2, labels = "s62 log2 FPKM", cex.axis = 1.5, at = 5, line = 1.5, tck = 0) #y

text(13.8, -5.1, labels=rho_S61_S62, cex = 1.2)
text(10,-5, "\\*r", vfont = c("sans serif", "bold italic"), cex=2)
text(11.5,-5, "=", vfont = c("sans serif", "plain"), cex = 2)
