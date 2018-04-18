#Boxplots A vs Z

setwd("G:/Masterarbete/DC_analysis")

instar_V <- read.delim("instar_V-assigned_A_or_Z-filtered.txt", header = TRUE)
pupa <- read.delim("pupa-assigned_A_or_Z-filtered.txt", header = TRUE)
adult <- read.delim("adult-assigned_A_or_Z-filtered.txt", header = TRUE)

instar_V_female <- instar_V[c(4,6)]
instar_V_female$group <- rep("1_instar_V_female", nrow(instar_V_female))
instar_V_female$group <- paste(instar_V_female$group, instar_V_female$chromosome)
instar_V_female <- instar_V_female[c(1,3)]
names(instar_V_female)[1] <-"FPKM"

pupa_female <- pupa[c(4,6)]
pupa_female$group <- rep("3_pupa_female", nrow(pupa_female))
pupa_female$group <- paste(pupa_female$group, pupa_female$chromosome)
pupa_female <- pupa_female[c(1,3)]
names(pupa_female)[1] <-"FPKM"

adult_female <- adult[c(4,6)]
adult_female$group <- rep("5_adult_female", nrow(adult_female))
adult_female$group <- paste(adult_female$group, adult_female$chromosome)
adult_female <- adult_female[c(1,3)]
names(adult_female)[1] <-"FPKM"

instar_V_male <- instar_V[c(5,6)]
instar_V_male$group <- rep("2_instar_V_male", nrow(instar_V_male))
instar_V_male$group <- paste(instar_V_male$group, instar_V_male$chromosome)
instar_V_male <- instar_V_male[c(1,3)]
names(instar_V_male)[1] <-"FPKM"

pupa_male <- pupa[c(5,6)]
pupa_male$group <- rep("4_pupa_male", nrow(pupa_male))
pupa_male$group <- paste(pupa_male$group, pupa_male$chromosome)
pupa_male <- pupa_male[c(1,3)]
names(pupa_male)[1] <-"FPKM"

adult_male <- adult[c(5,6)]
adult_male$group <- rep("6_adult_male", nrow(adult_male))
adult_male$group <- paste(adult_male$group, adult_male$chromosome)
adult_male <- adult_male[c(1,3)]
names(adult_male)[1] <-"FPKM"

all_samples <- rbind(instar_V_female, instar_V_male, pupa_female, pupa_male, adult_female, adult_male)

#Non-biased genes - sex-biased genes removed from filtered (FPKM>0) data set

instar_V_nonbiased <- read.delim("nonbiased_genes-instar_V-assigned_A_or_Z-filtered.txt", header = TRUE)
pupa_nonbiased <- read.delim("nonbiased_genes-pupa-assigned_A_or_Z-filtered.txt", header = TRUE)
adult_nonbiased <- read.delim("nonbiased_genes-adult-assigned_A_or_Z-filtered.txt", header = TRUE)

instar_V_female_nonbiased <- instar_V_nonbiased[c(4,6)]
instar_V_female_nonbiased$group <- rep("1_instar_V_female_nonbiased", nrow(instar_V_female_nonbiased))
instar_V_female_nonbiased$group <- paste(instar_V_female_nonbiased$group, instar_V_female_nonbiased$chromosome)
instar_V_female_nonbiased <- instar_V_female_nonbiased[c(1,3)]
names(instar_V_female_nonbiased)[1] <-"FPKM"

pupa_female_nonbiased <- pupa_nonbiased[c(4,6)]
pupa_female_nonbiased$group <- rep("3_pupa_female_nonbiased", nrow(pupa_female_nonbiased))
pupa_female_nonbiased$group <- paste(pupa_female_nonbiased$group, pupa_female_nonbiased$chromosome)
pupa_female_nonbiased <- pupa_female_nonbiased[c(1,3)]
names(pupa_female_nonbiased)[1] <-"FPKM"

adult_female_nonbiased <- adult_nonbiased[c(4,6)]
adult_female_nonbiased$group <- rep("5_adult_female_nonbiased", nrow(adult_female_nonbiased))
adult_female_nonbiased$group <- paste(adult_female_nonbiased$group, adult_female_nonbiased$chromosome)
adult_female_nonbiased <- adult_female_nonbiased[c(1,3)]
names(adult_female_nonbiased)[1] <-"FPKM"

instar_V_male_nonbiased <- instar_V_nonbiased[c(5,6)]
instar_V_male_nonbiased$group <- rep("2_instar_V_male_nonbiased", nrow(instar_V_male_nonbiased))
instar_V_male_nonbiased$group <- paste(instar_V_male_nonbiased$group, instar_V_male_nonbiased$chromosome)
instar_V_male_nonbiased <- instar_V_male_nonbiased[c(1,3)]
names(instar_V_male_nonbiased)[1] <-"FPKM"

pupa_male_nonbiased <- pupa_nonbiased[c(5,6)]
pupa_male_nonbiased$group <- rep("4_pupa_male_nonbiased", nrow(pupa_male_nonbiased))
pupa_male_nonbiased$group <- paste(pupa_male_nonbiased$group, pupa_male_nonbiased$chromosome)
pupa_male_nonbiased <- pupa_male_nonbiased[c(1,3)]
names(pupa_male_nonbiased)[1] <-"FPKM"

adult_male_nonbiased <- adult_nonbiased[c(5,6)]
adult_male_nonbiased$group <- rep("6_adult_male_nonbiased", nrow(adult_male_nonbiased))
adult_male_nonbiased$group <- paste(adult_male_nonbiased$group, adult_male_nonbiased$chromosome)
adult_male_nonbiased <- adult_male_nonbiased[c(1,3)]
names(adult_male_nonbiased)[1] <-"FPKM"

nonbiased <- rbind(instar_V_female_nonbiased, instar_V_male_nonbiased, pupa_female_nonbiased, pupa_male_nonbiased, adult_female_nonbiased, adult_male_nonbiased)


# Filtered FPKM > 0 ####

boxplot(log2(all_samples$FPKM)~all_samples$group, ylim = c(-8, 14),
        col = c("lightgrey", "darkorange"), notch = TRUE,
        at = c(1,2, 4,5,   8,9, 11,12,  15,16, 18,19),
        outline = FALSE, boxwex = 0.7,  xaxt = "n", frame.plot = FALSE)
axis(2, labels = "log2 FPKM(>0)", cex.axis = 1.2, at = 3, line = 1.5, tck = 0)

text(1.5, -8.5, "\\VE", vfont=c("sans serif","plain"), cex =2)
text(4.5, -8.5, "\\MA", vfont=c("sans serif","plain"), cex =2)

text(8.5, -8.5, "\\VE", vfont=c("sans serif","plain"), cex =2)
text(11.5, -8.5, "\\MA", vfont=c("sans serif","plain"), cex =2)

text(15.5, -8.5, "\\VE", vfont=c("sans serif","plain"), cex =2)
text(18.5, -8.5, "\\MA", vfont=c("sans serif","plain"), cex =2)

axis(1, labels = c("Instar V", "Pupa", "Adult"), lwd = 0, cex.axis = 1.5,
     at = c(3, 10, 17), line = 1)
legend(16.5, 15, legend = c("Autosomes","Z"), cex = 1,
       fill = c("lightgrey", "darkorange"), bty = "n")
segments(6.5, -7.5, 6.5, 10, lty = 3, lwd = 1)
segments(13.5, -7.5, 13.5, 10, lty = 3, lwd = 1)
segments(1, -7.5, 19, -7.5, lty = 1, lwd = 1)

##############################################################################

#Before vs after removing SBG####

boxplot(log2(all_samples$FPKM)~all_samples$group, ylim = c(-8, 14),
        col = c("grey90", "grey90"), notch = TRUE,
        at = c(1,2, 4,5,   8,9, 11,12,  15,16, 18,19),
        outline = FALSE, boxwex = 0.7,  xaxt = "n", frame.plot = FALSE, border="grey65",
        boxlty = 2)
axis(2, labels = "log2 FPKM(>0)", cex.axis = 1.2, at = 3, line = 1.5, tck = 0)

boxplot(log2(nonbiased$FPKM)~nonbiased$group, ylim = c(-8, 14), add = TRUE, 
        col = c("lightgrey", "darkorange"), notch = TRUE,
        at = c(1.2,2.2, 4.2,5.2,   8.2,9.2, 11.2,12.2,  15.2,16.2, 18.2,19.2),
        outline = FALSE, boxwex = 0.7,  xaxt = "n", frame.plot = FALSE)

text(1.5, -8.5, "\\VE", vfont=c("sans serif","plain"), cex =2)
text(4.5, -8.5, "\\MA", vfont=c("sans serif","plain"), cex =2)

text(8.5, -8.5, "\\VE", vfont=c("sans serif","plain"), cex =2)
text(11.5, -8.5, "\\MA", vfont=c("sans serif","plain"), cex =2)

text(15.5, -8.5, "\\VE", vfont=c("sans serif","plain"), cex =2)
text(18.5, -8.5, "\\MA", vfont=c("sans serif","plain"), cex =2)

axis(1, labels = c("Instar V", "Pupa", "Adult"), lwd = 0, cex.axis = 1.5,
     at = c(3, 10, 17), line = 1)
legend(16.5, 15, legend = c("Autosomes","Z"), cex = 1,
       fill = c("lightgrey", "darkorange"), bty = "n")
segments(6.5, -7.5, 6.5, 10, lty = 3, lwd = 1)
segments(13.5, -7.5, 13.5, 10, lty = 3, lwd = 1)
segments(1, -7.5, 19, -7.5, lty = 1, lwd = 1)
