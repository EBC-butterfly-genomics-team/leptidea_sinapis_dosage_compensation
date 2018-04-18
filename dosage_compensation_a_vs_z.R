#DC analysis - calculate ratios between A and Z median values

#Input: tables of gene expression as mean FPKM-values for each group ####

setwd("G:/Masterarbete/DC_analysis")

instar_V <- read.delim("instar_V-assigned_A_or_Z-filtered.txt", header = TRUE)
pupa <- read.delim("pupa-assigned_A_or_Z-filtered.txt", header = TRUE)
adult <- read.delim("adult-assigned_A_or_Z-filtered.txt", header = TRUE)

#instar_V <- read.delim("nonbiased_genes-instar_V-assigned_A_or_Z-filtered.txt", header = TRUE)
#pupa <- read.delim("nonbiased_genes-pupa-assigned_A_or_Z-filtered.txt", header = TRUE)
#adult <- read.delim("nonbiased_genes-adult-assigned_A_or_Z-filtered.txt", header = TRUE)

#instar_V <- read.delim("instar_V-assigned_A_or_Z.txt", header = TRUE)
#pupa <- read.delim("pupa-assigned_A_or_Z.txt", header = TRUE)
#adult <- read.delim("adult-assigned_A_or_Z.txt", header = TRUE)


instar_V_A <- rbind(instar_V[instar_V$chromosome == "A", ])
instar_V_Z <- rbind(instar_V[instar_V$chromosome == "Z", ])
pupa_A <- rbind(pupa[pupa$chromosome == "A", ])
pupa_Z <- rbind(pupa[pupa$chromosome == "Z", ])
adult_A <- rbind(adult[adult$chromosome == "A", ])
adult_Z <- rbind(adult[adult$chromosome == "Z", ])

########################################################################

# Median FPKM ####
instar_V_male_a_median <- median(instar_V_A$FPKM_instar_V_male)
instar_V_male_z_median <- median(instar_V_Z$FPKM_instar_V_male)
instar_V_female_a_median <- median(instar_V_A$FPKM_instar_V_female)
instar_V_female_z_median <- median(instar_V_Z$FPKM_instar_V_female)

pupa_male_a_median <- median(pupa_A$FPKM_pupa_male)
pupa_male_z_median <- median(pupa_Z$FPKM_pupa_male)
pupa_female_a_median <- median(pupa_A$FPKM_pupa_female)
pupa_female_z_median <- median(pupa_Z$FPKM_pupa_female)

adult_male_a_median <- median(adult_A$FPKM_adult_male)
adult_male_z_median <- median(adult_Z$FPKM_adult_male)
adult_female_a_median <- median(adult_A$FPKM_adult_female)
adult_female_z_median <- median(adult_Z$FPKM_adult_female)

# Mean of FPKM ####
instar_V_male_a_mean <- mean(instar_V_A$FPKM_instar_V_male)
instar_V_male_z_mean <- mean(instar_V_Z$FPKM_instar_V_male)
instar_V_female_a_mean <- mean(instar_V_A$FPKM_instar_V_female)
instar_V_female_z_mean <- mean(instar_V_Z$FPKM_instar_V_female)

pupa_male_a_mean <- mean(pupa_A$FPKM_pupa_male)
pupa_male_z_mean <- mean(pupa_Z$FPKM_pupa_male)
pupa_female_a_mean <- mean(pupa_A$FPKM_pupa_female)
pupa_female_z_mean <- mean(pupa_Z$FPKM_pupa_female)

adult_male_a_mean <- mean(adult_A$FPKM_adult_male)
adult_male_z_mean <- mean(adult_Z$FPKM_adult_male)
adult_female_a_mean <- mean(adult_A$FPKM_adult_female)
adult_female_z_mean <- mean(adult_Z$FPKM_adult_female)

#Table of median and mean FPKM

median_and_mean <- matrix(c(instar_V_female_a_median, pupa_female_a_median, adult_female_a_median,
                            instar_V_female_z_median, pupa_female_z_median, adult_female_z_median,
                            instar_V_male_a_median, pupa_male_a_median, adult_male_a_median,
                            instar_V_male_z_median, pupa_male_z_median, adult_male_z_median,
                            instar_V_female_a_mean, pupa_female_a_mean, adult_female_a_mean,
                            instar_V_female_z_mean, pupa_female_z_mean, adult_female_z_mean,
                            instar_V_male_a_mean, pupa_male_a_mean, adult_male_a_mean,
                            instar_V_male_z_mean, pupa_male_z_mean, adult_male_z_mean),
                          nrow = 3, ncol = 8)

rownames(median_and_mean) <- c("Instar V","Pupa","Adult")
colnames(median_and_mean) <- c("Female median A", "Female median Z", 
                               "Male median A", "Male median Z",
                               "Female mean A", "Female mean Z",
                               "Male mean A", "Male mean Z")

write.table(median_and_mean, file = "median_and_mean_FPKM.txt", sep = "\t", col.names = NA, row.names = TRUE)
#write.table(median_and_mean, file = "nonbiased_genes_median_and_mean_FPKM.txt", sep = "\t", col.names = NA, row.names = TRUE)
#write.table(median_and_mean, file = "unfiltered_median_and_mean_FPKM.txt", sep = "\t", col.names = NA, row.names = TRUE)


# Ratios of median expression - male Z:A, female Z:A, A m:f, Z m:f ####
instar_V_male_ratio <- instar_V_male_z_median/instar_V_male_a_median
instar_V_female_ratio <- instar_V_female_z_median/instar_V_female_a_median
pupa_male_ratio <- pupa_male_z_median/pupa_male_a_median
pupa_female_ratio <- pupa_female_z_median/pupa_female_a_median
adult_male_ratio <- adult_male_z_median/adult_male_a_median
adult_female_ratio <- adult_female_z_median/adult_female_a_median

instar_V_a_ratio <- instar_V_male_a_median/instar_V_female_a_median
instar_V_z_ratio <- instar_V_male_z_median/instar_V_female_z_median
pupa_a_ratio <- pupa_male_a_median/pupa_female_a_median
pupa_z_ratio <- pupa_male_z_median/pupa_female_z_median
adult_a_ratio <- adult_male_a_median/adult_female_a_median
adult_z_ratio <- adult_male_z_median/adult_female_z_median

ratios_of_median <- matrix(c(instar_V_female_ratio, pupa_female_ratio, adult_female_ratio,
                             instar_V_male_ratio, pupa_male_ratio, adult_male_ratio,
                             instar_V_a_ratio, pupa_a_ratio, adult_a_ratio,
                             instar_V_z_ratio, pupa_z_ratio, adult_z_ratio),
                           nrow = 3, ncol = 4)
rownames(ratios_of_median) <- c("Instar V", "Pupa", "Adult")
colnames(ratios_of_median) <- c("Female Z:A ratio", "Male Z:A ratio", "A M:F ratio", "Z M:F ratio")

write.table(ratios_of_median, file = "ratios_of_median_FPKM.txt", sep = "\t", col.names = NA, row.names = TRUE)
#write.table(ratios_of_median, file = "nonbiased_genes_ratios_of_median_FPKM.txt", sep = "\t", col.names = NA, row.names = TRUE)


# MWU-test of significant difference between Z and A for each sex, and Z M/F, A M/F ####
instar_V_male_mwu <- wilcox.test(instar_V_A$FPKM_instar_V_male, instar_V_Z$FPKM_instar_V_male,
                                 paired = FALSE)
instar_V_female_mwu <- wilcox.test(instar_V_A$FPKM_instar_V_female, instar_V_Z$FPKM_instar_V_female,
                                 paired = FALSE)
pupa_male_mwu <- wilcox.test(pupa_A$FPKM_pupa_male, pupa_Z$FPKM_pupa_male,
                                 paired = FALSE)
pupa_female_mwu <- wilcox.test(pupa_A$FPKM_pupa_female, pupa_Z$FPKM_pupa_female,
                             paired = FALSE)
adult_male_mwu <- wilcox.test(adult_A$FPKM_adult_male, adult_Z$FPKM_adult_male,
                             paired = FALSE)
adult_female_mwu <- wilcox.test(adult_A$FPKM_adult_female, adult_Z$FPKM_adult_female,
                              paired = FALSE)

instar_V_a_mwu <-wilcox.test(instar_V_A$FPKM_instar_V_male, instar_V_A$FPKM_instar_V_female,
                             paired = FALSE)
instar_V_z_mwu <-wilcox.test(instar_V_Z$FPKM_instar_V_male, instar_V_Z$FPKM_instar_V_female,
                             paired = FALSE)
pupa_a_mwu <-wilcox.test(pupa_A$FPKM_pupa_male, pupa_A$FPKM_pupa_female,
                             paired = FALSE)
pupa_z_mwu <-wilcox.test(pupa_Z$FPKM_pupa_male, pupa_Z$FPKM_pupa_female,
                             paired = FALSE)
adult_a_mwu <-wilcox.test(adult_A$FPKM_adult_male, adult_A$FPKM_adult_female,
                         paired = FALSE)
adult_z_mwu <-wilcox.test(adult_Z$FPKM_adult_male, adult_Z$FPKM_adult_female,
                         paired = FALSE)

mwu <- matrix(c(instar_V_female_mwu$p.value, pupa_female_mwu$p.value, adult_female_mwu$p.value,
                instar_V_male_mwu$p.value, pupa_male_mwu$p.value, adult_male_mwu$p.value,
                instar_V_a_mwu$p.value, pupa_a_mwu$p.value, adult_a_mwu$p.value,
                instar_V_z_mwu$p.value, pupa_z_mwu$p.value, adult_z_mwu$p.value),
              nrow = 3, ncol = 4)
rownames(mwu) <- c("Instar V", "Pupa", "Adult")
colnames(mwu) <- c("Female - A vs Z", "Male - A vs Z", "A - F vs M", "Z - F vs M")

write.table(mwu, file = "mwu_test_p-values.txt", sep = "\t", col.names = NA, row.names = TRUE)
#write.table(mwu, file = "nonbiased_genes_mwu_test_p-values.txt", sep = "\t", col.names = NA, row.names = TRUE)


# Bootstraping for median ratio confidence intervals ####

library(boot)

# Adult Autosome M:F ratio ####
adult_a_median_ratio_function <- function(data, indices) {
  adult_A <- data[indices,]
  (median(adult_A$FPKM_adult_male))/(median(adult_A$FPKM_adult_female))
}

adult_a_median_boot.res <- boot(data = adult_A, statistic = adult_a_median_ratio_function, R = 10000)
print(adult_a_median_boot.res)
plot(adult_a_median_boot.res)

adult_a_median_CI <- boot.ci(boot.out = adult_a_median_boot.res, type = "basic")
print(adult_a_median_CI)

# Adult Z-chromosome M:F ratio ####

adult_z_median_ratio_function <- function(data, indices) {
  adult_Z <- data[indices,]
  (median(adult_Z$FPKM_adult_male))/(median(adult_Z$FPKM_adult_female))
}

adult_z_median_boot.res <- boot(data = adult_Z, statistic = adult_z_median_ratio_function, R = 10000)
print(adult_z_median_boot.res)
plot(adult_z_median_boot.res)

adult_z_median_CI <- boot.ci(boot.out = adult_z_median_boot.res, type = "basic")
print(adult_z_median_CI)

# Pupa Autosome M:F ratio ####
pupa_a_median_ratio_function <- function(data, indices) {
  pupa_A <- data[indices,]
  (median(pupa_A$FPKM_pupa_male))/(median(pupa_A$FPKM_pupa_female))
}

pupa_a_median_boot.res <- boot(data = pupa_A, statistic = pupa_a_median_ratio_function, R = 10000)
print(pupa_a_median_boot.res)
plot(pupa_a_median_boot.res)

pupa_a_median_CI <- boot.ci(boot.out = pupa_a_median_boot.res, type = "basic")
print(pupa_a_median_CI)

# Pupa Z-chromosome M:F ratio ####
pupa_z_median_ratio_function <- function(data, indices) {
  pupa_Z <- data[indices,]
  (median(pupa_Z$FPKM_pupa_male))/(median(pupa_Z$FPKM_pupa_female))
}

pupa_z_median_boot.res <- boot(data = pupa_Z, statistic = pupa_z_median_ratio_function, R = 10000)
print(pupa_z_median_boot.res)
plot(pupa_z_median_boot.res)

pupa_z_median_CI <- boot.ci(boot.out = pupa_z_median_boot.res, type = "basic")
print(pupa_z_median_CI)


# Instar V Autosome M:F ratio ####
instar_V_a_median_ratio <- function(data, indices) {
  instar_V_A <- data[indices,]
  (median(instar_V_A$FPKM_instar_V_male))/(median(instar_V_A$FPKM_instar_V_female))
}

instar_V_a_median_boot.res <- boot(data = instar_V_A, statistic = instar_V_a_median_ratio, R = 10000)
print(instar_V_a_median_boot.res)
plot(instar_V_a_median_boot.res)

instar_V_a_median_CI <- boot.ci(boot.out = instar_V_a_median_boot.res, type = "basic")
print(instar_V_a_median_CI)

# Instar V Z-chromosome M:F ratio ####
instar_V_z_median_ratio <- function(data, indices) {
  instar_V_Z <- data[indices,]
  (median(instar_V_Z$FPKM_instar_V_male))/(median(instar_V_Z$FPKM_instar_V_female))
}

instar_V_z_median_boot.res <- boot(data = instar_V_Z, statistic = instar_V_z_median_ratio, R = 10000)
print(instar_V_z_median_boot.res)
plot(instar_V_z_median_boot.res)

instar_V_z_median_CI <- boot.ci(boot.out = instar_V_z_median_boot.res, type = "basic")
print(instar_V_z_median_CI)

########################################

# Adult Male Z:A ratio ####
adult_male_median_function <- function(data, indices) {
  adult <- data[indices,]
  (median(adult$FPKM_adult_male[adult$chromosome=="Z"]))/(median(adult$FPKM_adult_male[adult$chromosome=="A"]))
}

adult_male_median_boot <- boot(data = adult, statistic = adult_male_median_function, R = 10000)
print(adult_male_median_boot)
plot(adult_male_median_boot)

adult_male_median_CI <- boot.ci(boot.out = adult_male_median_boot, type = "basic")
print(adult_male_median_CI)

# Adult Female Z:A ratio ####
adult_female_median_function <- function(data, indices) {
  adult <- data[indices,]
  (median(adult$FPKM_adult_female[adult$chromosome=="Z"]))/(median(adult$FPKM_adult_female[adult$chromosome=="A"]))
}

adult_female_median_boot <- boot(data = adult, statistic = adult_female_median_function, R = 10000)
print(adult_female_median_boot)
plot(adult_female_median_boot)

adult_female_median_CI <- boot.ci(boot.out = adult_female_median_boot, type = "basic")
print(adult_female_median_CI)

# Pupa Male Z:A ratio ####
pupa_male_median_function <- function(data, indices) {
  pupa <- data[indices,]
  (median(pupa$FPKM_pupa_male[pupa$chromosome=="Z"]))/(median(pupa$FPKM_pupa_male[pupa$chromosome=="A"]))
}

pupa_male_median_boot <- boot(data = pupa, statistic = pupa_male_median_function, R = 10000)
print(pupa_male_median_boot)
plot(pupa_male_median_boot)

pupa_male_median_CI <- boot.ci(boot.out = pupa_male_median_boot, type = "basic")
print(pupa_male_median_CI)

# Pupa Female Z:A ratio ####
pupa_female_median_function <- function(data, indices) {
  pupa <- data[indices,]
  (median(pupa$FPKM_pupa_female[pupa$chromosome=="Z"]))/(median(pupa$FPKM_pupa_female[pupa$chromosome=="A"]))
}

pupa_female_median_boot <- boot(data = pupa, statistic = pupa_female_median_function, R = 10000)
print(pupa_female_median_boot)
plot(pupa_female_median_boot)

pupa_female_median_CI <- boot.ci(boot.out = pupa_female_median_boot, type = "basic")
print(pupa_female_median_CI)


# Instar V Male Z:A ratio ####
instar_V_male_median_function <- function(data, indices) {
  instar_V <- data[indices,]
  (median(instar_V$FPKM_instar_V_male[instar_V$chromosome=="Z"]))/(median(instar_V$FPKM_instar_V_male[instar_V$chromosome=="A"]))
}

instar_V_male_median_boot <- boot(data = instar_V, statistic = instar_V_male_median_function, R = 10000)
print(instar_V_male_median_boot)
plot(instar_V_male_median_boot)

instar_V_male_median_CI <- boot.ci(boot.out = instar_V_male_median_boot, type = "basic")
print(instar_V_male_median_CI)

# Instar V Female Z:A ratio ####
instar_V_female_median_function <- function(data, indices) {
  instar_V <- data[indices,]
  (median(instar_V$FPKM_instar_V_female[instar_V$chromosome=="Z"]))/(median(instar_V$FPKM_instar_V_female[instar_V$chromosome=="A"]))
}

instar_V_female_median_boot <- boot(data = instar_V, statistic = instar_V_female_median_function, R = 10000)
print(instar_V_female_median_boot)
plot(instar_V_female_median_boot)

instar_V_female_median_CI <- boot.ci(boot.out = instar_V_female_median_boot, type = "basic")
print(instar_V_female_median_CI)

##############################################