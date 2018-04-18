#!/bin/bash

#Lars Höök 2017

#Sex-biased genes are removed from genes filtered FPKM>0 - for dosage compensation analysis
#For Fisher test (expected distribution), ALL genes filtered in DEseq2 are assign as A or Z (including SBG)
#Male and female biased genes are picked from list of all genes assigned A or Z, but not filtered FPKM>0
#Gene names are added to sex biased genes for volcano plot

STAGE1="instar_V-assigned_A_or_Z.txt"
FILTERED1="instar_V-assigned_A_or_Z-filtered.txt"
MBG1="MBG_Instar.txt"
FBG1="FBG_Instar.txt"
SBG_STAGE1="Filtered_LFC_1.5_baseMean_10_Instar.txt"
ALL_DESEQ1="LFC_Instar.txt"
GENE_NAMES1="instar_V-concatenated.txt"

STAGE2="pupa-assigned_A_or_Z.txt"
FILTERED2="pupa-assigned_A_or_Z-filtered.txt"
MBG2="MBG_Pupa.txt"
FBG2="FBG_Pupa.txt"
SBG_STAGE2="Filtered_LFC_1.5_baseMean_10_Pupa.txt"
ALL_DESEQ2="LFC_Pupa.txt"
GENE_NAMES2="pupa-concatenated.txt"

STAGE3="adult-assigned_A_or_Z.txt"
FILTERED3="adult-assigned_A_or_Z-filtered.txt"
MBG3="MBG_Adult.txt"
FBG3="FBG_Adult.txt"
SBG_STAGE3="Filtered_LFC_1.5_baseMean_10_Adult.txt"
ALL_DESEQ3="LFC_Adult.txt"
GENE_NAMES3="adult-concatenated.txt"


sort $MBG1 -o $MBG1
sort $FBG1 -o $FBG1

grep -vf <(cut -f1 $SBG_STAGE1 | sort) $FILTERED1 > nonbiased_genes-$FILTERED1
grep -f <(cut -f1 $ALL_DESEQ1 | sort) $STAGE1 > FisherTest-$STAGE1
grep -f <(cut -f1 $MBG1) $STAGE1 > FisherTest_MBG-$STAGE1
grep -f <(cut -f1 $FBG1) $STAGE1 > FisherTest_FBG-$STAGE1

paste <(grep -f <(cut -f1 $MBG1) <(cut -f1,2 $GENE_NAMES1)) <(cut -f2-8 $MBG1) > gene_names-$MBG1
paste <(grep -f <(cut -f1 $FBG1) <(cut -f1,2 $GENE_NAMES1)) <(cut -f2-8 $FBG1) > gene_names-$FBG1


sort $MBG2 -o $MBG2
sort $FBG2 -o $FBG2

grep -vf <(cut -f1 $SBG_STAGE2 | sort) $FILTERED2 > nonbiased_genes-$FILTERED2
grep -f <(cut -f1 $ALL_DESEQ2 | sort) $STAGE2 > FisherTest-$STAGE2
grep -f <(cut -f1 $MBG2 | sort) $STAGE2 > FisherTest_MBG-$STAGE2
grep -f <(cut -f1 $FBG2 | sort) $STAGE2 > FisherTest_FBG-$STAGE2

paste <(grep -f <(cut -f1 $MBG2) <(cut -f1,2 $GENE_NAMES2)) <(cut -f2-8 $MBG2) > gene_names-$MBG2
paste <(grep -f <(cut -f1 $FBG2) <(cut -f1,2 $GENE_NAMES2)) <(cut -f2-8 $FBG2) > gene_names-$FBG2


sort $MBG3 -o $MBG3
sort $FBG3 -o $FBG3

grep -vf <(cut -f1 $SBG_STAGE3 | sort) $FILTERED3 > nonbiased_genes-$FILTERED3
grep -f <(cut -f1 $ALL_DESEQ3 | sort) $STAGE3 > FisherTest-$STAGE3
grep -f <(cut -f1 $MBG3 | sort) $STAGE3 > FisherTest_MBG-$STAGE3
grep -f <(cut -f1 $FBG3 | sort) $STAGE3 > FisherTest_FBG-$STAGE3

paste <(grep -f <(cut -f1 $MBG3) <(cut -f1,2 $GENE_NAMES3)) <(cut -f2-8 $MBG3) > gene_names-$MBG3
paste <(grep -f <(cut -f1 $FBG3) <(cut -f1,2 $GENE_NAMES3)) <(cut -f2-8 $FBG3) > gene_names-$FBG3

