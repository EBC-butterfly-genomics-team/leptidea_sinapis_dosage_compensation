#Boxplots of individual chromosomes

setwd("G:/Masterarbete/DC_analysis")

instar_V <- read.delim("instar_V-assigned_to_chromosomes-filtered.txt", header = TRUE)
pupa <- read.delim("pupa-assigned_to_chromosomes-filtered.txt", header = TRUE)
adult <- read.delim("adult-assigned_to_chromosomes-filtered.txt", header = TRUE)

instar_V_female <- instar_V[c(4,6)]
instar_V_female_z <- rbind(instar_V_female[instar_V_female$chromosome=="0", ])
instar_V_female_z_median <- median(instar_V_female_z$FPKM_instar_V_female)

instar_V_male <- instar_V[c(5,6)]
instar_V_male_z <- rbind(instar_V_male[instar_V_male$chromosome=="0", ])
instar_V_male_z_median <- median(instar_V_male_z$FPKM_instar_V_male)

pupa_female <- pupa[c(4,6)]
pupa_female_z <- rbind(pupa_female[pupa_female$chromosome=="0", ])
pupa_female_z_median <- median(pupa_female_z$FPKM_pupa_female)

pupa_male <- pupa[c(5,6)]
pupa_male_z <- rbind(pupa_male[pupa_male$chromosome=="0", ])
pupa_male_z_median <- median(pupa_male_z$FPKM_pupa_male)

adult_female <- adult[c(4,6)]
adult_female_z <- rbind(adult_female[adult_female$chromosome=="0", ])
adult_female_z_median <- median(adult_female_z$FPKM_adult_female)

adult_male <- adult[c(5,6)]
adult_male_z <- rbind(adult_male[adult_male$chromosome=="0", ])
adult_male_z_median <- median(adult_male_z$FPKM_adult_male)



####################################################################

par(pty="s", mfrow = c(1,3))

boxplot(log2(instar_V_female$FPKM_instar_V_female)~instar_V_female$chromosome, ylim = c(-8, 14), notch = TRUE,
        col = c(rep("darkorange",1),rep("grey",20)), 
        outline = FALSE, xaxt = "n", boxwex = 0.7, frame.plot = FALSE,
        main = "Instar V", cex.main = 2)
axis(3, labels = "A", cex.axis = 3, at = 0, tck = 0)
axis(2, labels = "log2 FPKM(>0)", cex.axis = 1.2, at = 3, line = 1.5, tck = 0)
axis(1, las=2, cex.axis = 1.3, labels = c("Z", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"),
     at = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21))

segments(1, log2(instar_V_female_z_median), 21.5, log2(instar_V_female_z_median), col = "darkorange", lty = 2, lwd = 2)

boxplot(log2(pupa_female$FPKM_pupa_female)~pupa_female$chromosome, ylim = c(-8, 14), notch = TRUE,
        col = c(rep("darkorange",1),rep("grey",20)), 
        outline = FALSE, xaxt = "n", boxwex = 0.7, frame.plot = FALSE,
        main = "Pupa", cex.main = 2)
axis(1, las=2, cex.axis = 1.3, labels = c("Z", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"),
     at = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21))

segments(1, log2(pupa_female_z_median), 21.5, log2(pupa_female_z_median), col = "darkorange", lty = 2, lwd = 2)

axis(1, labels = c("Female chromosomes"), lwd = 0, at = 10.5, cex.axis = 1.5, tck = 0, line = 3)

boxplot(log2(adult_female$FPKM_adult_female)~adult_female$chromosome, ylim = c(-8, 14), notch = TRUE,
        col = c(rep("darkorange",1),rep("grey",20)), 
        outline = FALSE, xaxt = "n", boxwex = 0.7, frame.plot = FALSE,
        main = "Adult", cex.main = 2)
axis(1, las=2, cex.axis = 1.3, labels = c("Z", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"),
     at = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21))

segments(1, log2(adult_female_z_median), 21.5, log2(adult_female_z_median), col = "darkorange", lty = 2, lwd = 2)

#############################################################################

par(pty="s", mfrow = c(1,3))

boxplot(log2(instar_V_male$FPKM_instar_V_male)~instar_V_male$chromosome, ylim = c(-8, 14), notch = TRUE,
        col = c(rep("darkorange",1),rep("grey",20)), 
        outline = FALSE, xaxt = "n", boxwex = 0.7, frame.plot = FALSE,
        main = "Instar V", cex.main = 2)
axis(2, labels = "log2 FPKM(>0)", cex.axis = 1.2, at = 3, line = 1.5, tck = 0)
axis(3, labels = "B", cex.axis = 3, at = 0, tck = 0)
axis(1, las = 2, cex.axis = 1.3, labels = c("Z", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"),
     at = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21))

segments(1, log2(instar_V_male_z_median), 21.5, log2(instar_V_male_z_median), col = "darkorange", lty = 2, lwd = 2)

boxplot(log2(pupa_male$FPKM_pupa_male)~pupa_male$chromosome, ylim = c(-8, 14), notch = TRUE,
        col = c(rep("darkorange",1),rep("grey",20)), 
        outline = FALSE, xaxt = "n", boxwex = 0.7, frame.plot = FALSE,
        main = "Pupa", cex.main = 2)
axis(1, las=2, cex.axis = 1.3, labels = c("Z", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"),
     at = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21))

segments(1, log2(pupa_male_z_median), 21.5, log2(pupa_male_z_median), col = "darkorange", lty = 2, lwd = 2)

axis(1, labels = c("Male chromosomes"), lwd = 0, at = 10.5, cex.axis = 1.5, tck = 0, line = 3)

boxplot(log2(adult_male$FPKM_adult_male)~adult_male$chromosome, ylim = c(-8, 14), notch = TRUE,
        col = c(rep("darkorange",1),rep("grey",20)), 
        outline = FALSE, xaxt = "n", boxwex = 0.7, frame.plot = FALSE,
        main = "Adult", cex.main = 2)
axis(1, las=2, cex.axis = 1.3, labels = c("Z", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"),
     at = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21))

segments(1, log2(adult_male_z_median), 21.5, log2(adult_male_z_median), col = "darkorange", lty = 2, lwd = 2)
