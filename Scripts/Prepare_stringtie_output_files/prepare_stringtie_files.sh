#!/bin/bash


##################################################################

################################################

# 1: Path to StringTie main folder and filename:

################################################

StringTie=stringtie				#Main folder containing StringTie output subfolders for each sample
STAR=STAR.Aligned.out.gene_abund.tab		#STAR-StringTie outfile name (same for all samples)


#######################################

# 2: List with assigned scaffolds

#######################################

ASSIGNED_A_OR_Z=mapped_scaffolds_filtered_above_200_assigned_a_or_z_by_0.95.txt
ASSIGNED_CHR=mapped_scaffolds_filtered_above_200_assigned_to_chromosome_by_0.9.txt


#############################

# 3: Name groups and samples:

#############################

#Groups are combined two-by-two for filtering and comparison (G1 vs G2 etc.)

STAGE1="instar_V"

GROUP1="instar_V_female"
S1=S49
S2=S52
S3=S56

GROUP2="instar_V_male"
S4=S48 
S5=S51
S6=S55

STAGE2="pupa"

GROUP3="pupa_female"
S7=S35
S8=S53
S9=S58

GROUP4="pupa_male"
S10=S34
S11=S36
S12=S57

STAGE3="adult"

GROUP5="adult_female"
S13=S59
S14=S38
S15=S10

GROUP6="adult_male"
S16=S60
S17=S61
S18=S62


############################################

# 4:Run script: ./prepare_stringtie_files.sh

############################################


date > prepare_stringtie_files.log


########################

mkdir -p samples/$GROUP1

echo
echo preparing $GROUP1 samples...

perl 1_sort_samples.pl $StringTie/*$S1/$STAR samples/$GROUP1/$S1.gene_abund.tab
perl 1_sort_samples.pl $StringTie/*$S2/$STAR samples/$GROUP1/$S2.gene_abund.tab
perl 1_sort_samples.pl $StringTie/*$S3/$STAR samples/$GROUP1/$S3.gene_abund.tab

echo calulating mean FPKM...

cat samples/$GROUP1/*.tab | perl 2_calculate_mean.pl samples/$GROUP1/FPKM_$GROUP1.txt $GROUP1

mkdir -p samples/$GROUP2

echo preparing $GROUP2 samples...

perl 1_sort_samples.pl $StringTie/*$S4/$STAR samples/$GROUP2/$S4.gene_abund.tab
perl 1_sort_samples.pl $StringTie/*$S5/$STAR samples/$GROUP2/$S5.gene_abund.tab
perl 1_sort_samples.pl $StringTie/*$S6/$STAR samples/$GROUP2/$S6.gene_abund.tab

echo calulating mean FPKM...

cat samples/$GROUP2/*.tab | perl 2_calculate_mean.pl samples/$GROUP2/FPKM_$GROUP2.txt $GROUP2

echo concatenating $STAGE1...
paste samples/$GROUP1/FPKM_$GROUP1.txt samples/$GROUP2/FPKM_$GROUP2.txt | cut -f 1-4,8 > $STAGE1-concatenated.txt

echo assigning genes as A or Z...
perl 3a_assign_genes.pl $ASSIGNED_A_OR_Z $STAGE1-concatenated.txt

echo assigning genes to chromosomes...
perl 3b_assign_genes_to_chromosomes.pl $ASSIGNED_CHR $STAGE1-concatenated.txt 

echo filtering out non-expressed genes...
perl 4_filter_genes.pl $STAGE1-assigned_A_or_Z.txt assigned_A_or_Z $STAGE1-assigned_A_or_Z-filtered.txt
perl 4_filter_genes.pl $STAGE1-assigned_to_chromosomes.txt assigned_to_chromosomes $STAGE1-assigned_to_chromosomes-filtered.txt


########################


mkdir -p samples/$GROUP3

echo
echo preparing $GROUP3 samples...

perl 1_sort_samples.pl $StringTie/*$S7/$STAR samples/$GROUP3/$S7.gene_abund.tab
perl 1_sort_samples.pl $StringTie/*$S8/$STAR samples/$GROUP3/$S8.gene_abund.tab
perl 1_sort_samples.pl $StringTie/*$S9/$STAR samples/$GROUP3/$S9.gene_abund.tab

echo calculating mean FPKM...

cat samples/$GROUP3/*.tab | perl 2_calculate_mean.pl samples/$GROUP3/FPKM_$GROUP3.txt $GROUP3

mkdir -p samples/$GROUP4

echo preparing $GROUP4 samples...

perl 1_sort_samples.pl $StringTie/*$S10/$STAR samples/$GROUP4/$S10.gene_abund.tab
perl 1_sort_samples.pl $StringTie/*$S11/$STAR samples/$GROUP4/$S11.gene_abund.tab
perl 1_sort_samples.pl $StringTie/*$S12/$STAR samples/$GROUP4/$S12.gene_abund.tab

echo calculating mean FPKM...

cat samples/$GROUP4/*.tab | perl 2_calculate_mean.pl samples/$GROUP4/FPKM_$GROUP4.txt $GROUP4

echo concatenating $STAGE2...
paste samples/$GROUP3/FPKM_$GROUP3.txt samples/$GROUP4/FPKM_$GROUP4.txt | cut -f 1-4,8 > $STAGE2-concatenated.txt

echo assigning genes as A or Z...
perl 3a_assign_genes.pl $ASSIGNED_A_OR_Z $STAGE2-concatenated.txt

echo assigning genes to chromosomes...
perl 3b_assign_genes_to_chromosomes.pl $ASSIGNED_CHR $STAGE2-concatenated.txt 

echo filtering out non-expressed genes...
perl 4_filter_genes.pl $STAGE2-assigned_A_or_Z.txt assigned_A_or_Z $STAGE2-assigned_A_or_Z-filtered.txt
perl 4_filter_genes.pl $STAGE2-assigned_to_chromosomes.txt assigned_to_chromosomes $STAGE2-assigned_to_chromosomes-filtered.txt


########################


mkdir -p samples/$GROUP5

echo
echo preparing $GROUP5 samples...

perl 1_sort_samples.pl $StringTie/*$S13/$STAR samples/$GROUP5/$S13.gene_abund.tab
perl 1_sort_samples.pl $StringTie/*$S14/$STAR samples/$GROUP5/$S14.gene_abund.tab
perl 1_sort_samples.pl $StringTie/*$S15/$STAR samples/$GROUP5/$S15.gene_abund.tab

echo calculating mean FPKM...
cat samples/$GROUP5/*.tab | perl 2_calculate_mean.pl samples/$GROUP5/FPKM_$GROUP5.txt $GROUP5

mkdir -p samples/$GROUP6

echo preparing $GROUP6 samples...

perl 1_sort_samples.pl $StringTie/*$S16/$STAR samples/$GROUP6/$S16.gene_abund.tab
perl 1_sort_samples.pl $StringTie/*$S17/$STAR samples/$GROUP6/$S17.gene_abund.tab
perl 1_sort_samples.pl $StringTie/*$S18/$STAR samples/$GROUP6/$S18.gene_abund.tab

echo calculating mean FPKM...
cat samples/$GROUP6/*.tab | perl 2_calculate_mean.pl samples/$GROUP6/FPKM_$GROUP6.txt $GROUP6

echo concatenating $STAGE3...
paste samples/$GROUP5/FPKM_$GROUP5.txt samples/$GROUP6/FPKM_$GROUP6.txt | cut -f 1-4,8 > $STAGE3-concatenated.txt

echo assigning genes as A or Z...
perl 3a_assign_genes.pl $ASSIGNED_A_OR_Z $STAGE3-concatenated.txt

echo assigning genes to chromosomes...
perl 3b_assign_genes_to_chromosomes.pl $ASSIGNED_CHR $STAGE3-concatenated.txt 

echo filtering out non-expressed genes...
perl 4_filter_genes.pl $STAGE3-assigned_A_or_Z.txt assigned_A_or_Z $STAGE3-assigned_A_or_Z-filtered.txt
perl 4_filter_genes.pl $STAGE3-assigned_to_chromosomes.txt assigned_to_chromosomes $STAGE3-assigned_to_chromosomes-filtered.txt



#############################################################################
