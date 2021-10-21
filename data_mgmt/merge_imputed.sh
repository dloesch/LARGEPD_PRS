#!/usr/bin/bash
#$ -q threaded.q
#$ -P toconnor-lab
#$ -pe thread 4
#$ -l mem_free=10G
#$ -e merge.log
#$ -o merge.log
#$ -cwd
#$ -N merge


###plink
export PATH=/usr/local/packages/plink-1.90.beta-3.6/bin:$PATH

#specify chrom
chr=$1

tb=/local/chib/oconnor_genomes/GLAD/data_after_imputation/luoetal_TBprogression/chr$chr.dose.vcf.gz
large=/local/chib/toconnor_grp/LARGE-PD/imputed_2021/chr$chr.dose.vcf.gz

#merge
bcftools merge -m none -Oz -o tb_large.chr$chr.vcf.gz --threads 4 $tb $large

#ld prune
plink --vcf tb_large.chr$chr.vcf.gz --indep-pairwise 50 5 0.2 --out ld.chr$chr --maf 0.01 --geno 0.0

#extract
plink --vcf tb_large.chr$chr.vcf.gz --extract ld.chr$chr.prune.in --out tb_large.chr$chr.pruned --recode vcf-iid bgz --output-chr chr26
