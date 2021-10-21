#!/usr/bin/bash
#$ -q all.q  
#$ -P toconnor-lab
#$ -l mem_free=5G
#$ -e PRS.log
#$ -o PRS.log
#$ -cwd
#$ -N PRS



#run PRS scripts


###calculate PRS
sumstats=nalls.PD.sumstats.txt

#genotype file
#geno=large-PD.files.txt
geno=1KG.files.txt
#geno=tb.files.txt
#geno=pgp.files.txt

#prefix
#prefix=large
prefix=1KG
#prefix=tb
#prefix=pgp

Rscript ./scripts/calc_PRS.R -s $sumstats -g $geno -p $prefix -b HG38

##test
pheno=large.pheno.03_2021.txt
PRS=$prefix.PRS.txt
trait=PD_STATUS
covars=AGE,SEX,PC1,PC2,PC3,P4,PC5,SITE
k=0.005
#Rscript test_PRS.R --PRS $PRS --pheno $pheno -t $trait --prevalence $k -c $covars
