module load PRSice
module load plink
module load R


directory=./Doug/PRS
base=$directory/PD.sumstats.txt
target=$directory/large.PRS #plink fileset
pheno=$directory/large.pheno.txt
ld_ref=$directory/1KGP.EUR #plink fileset of 1KGP Europeans with a MAF of 0.05

#if using the provided Rscript wrapper: 
#Rscript $directory/PRSice.R --dir $directory --prsice PRSice
 PRSice \
    --base $base \
    --snp SNP --chr CHR --bp BP \
    --A1 A1 --A2 A2 \
    --pvalue P --stat BETA \
    --target $target \
    --pheno $pheno --ignore-fid --pheno-col DISEASE \
    --ld  $ld_ref \
    --binary-target T \
    --perm 10000 \
     --thread 2 \
    --prevalence 0.005 \
    --out 1KGP.$pop.full \
    --cov $pheno --cov-col AGE,SEX,SITE,@[PC1-10] --cov-factor SEX,SITE \
    