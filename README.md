# LARGEPD_PRS
 Scripts for the LARGE-PD PRS manuscript

## PRS
 Scripts for calculating and evaluating PRS models
1. calc_PRS.R: script for calculating the GWAS-significant PRS
2. PRS_functions.R: contains functions for calculating a PRS
3. test_PRS.R: script for testing a PRS using logistic regression and 10 folds
4. survival_analysis.R: script for generating K-M curves and Cox regression
5. PRS_plots.R: script for generating AUC plots and OR plots
6. plot_PRS_dist.R: script for plotting distribution of PD PRS
7. PRS_wilcox.R: script for characterizing the PD PRS by LARGE-PD recruitment site
8. PD_PRS_1KG.R: script for charaterizing the PD PRS in 1000 Genomes
9. AUC_PVAL.R: script for getting pvalues for the AUC of different PRS models
10. PRS_PCA_plots.R: script for plotting PCA and PRS distributions
11. PRS.sh: sample bash script for calculating and testing PRS
12. run_PRSice.sh: bash script for running PRSice-2

## haplotypes
 Scripts for performing haplotype analysis in SNCA
1. hap_blocks.sh: bash script for using PLINK to estimate haplotype blocks
2. haplotpes.R: script for extracting haplotypes from a phased VCF file
3. hap_functions.R: functions for extracting haplotypes
4. plot_haplotypes.R: script for plotting haplotypes
5. test_haplotypes.R: script for testing haplotyeps for association with PD
6. haplotype_descriptive_stats.R: script for getting descriptive stats from haplotype data
7. convert_nexus.R: script for converting haplotpye data to nexus format for visualizing in PopArt

## data_mgmt
 Data management scripts
1. align_merge.sh: script for aligning and merging all datasets used in study
2. subset_SNCA.sh: script for subsetting region around rs356182 in SNCA for datasets used in study
3. merge_imputed.sh: script for merging imputed data from LARGE-PD and Luo et al.
4. concat.sh: script for concatennating LARGE-PD and Luo et al. data
5. align: contains scipts for aligning genotype data, with align.R being the primary script. 

## misc
 miscellaneous scripts used in this project
1. filter_rels.R: script for resolving relative pairs
2. get_PCSs.R: script for calculating PCs and a kinship matrix using PC-AiR and PC-Relate
3. merged_GWAS.R: script for running a GWAS usng GENESIS on the merged LARGE-PD/Luo et al. data
4. phase_beagle.sh : script for phasing using BEAGLE