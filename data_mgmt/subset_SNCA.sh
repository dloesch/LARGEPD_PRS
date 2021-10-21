#!/usr/bin/bash
#$ -q threaded.q  
#$ -P toconnor-lab
#$ -l mem_free=5G
#$ -e subset.log
#$ -o subset.log
#$ -pe thread 4
#$ -cwd
#$ -N subset


##script subsets out SNCA region based on file provided by IPDGC

start=$(grep -v "#"  data/NACHO_SNCA_subset.pvar | cut -f2 | head -1)
end=$(grep -v "#"  data/NACHO_SNCA_subset.pvar | cut -f2 | tail -1)

#1000 genomes
#IDs aren't set in this file
kg=/local/chib/toconnor_grp/1000G_GRC38_HighCov/CCDG_13607_B01_GRM_WGS_2019-02-19_chr4.recalibrated_variants.vcf.gz
bcftools view -r chr4:${start}-${end} -m2 -M2 -v snps -Oz -o temp.SNCA.vcf.gz --threads 4 $kg

bcftools annotate --set-id +'%CHROM\:%POS\:%REF\:%FIRST_ALT' -Oz -o 1KG.SNCA.vcf.gz --threads 4 temp.SNCA.vcf.gz

rm temp.SNCA.vcf.gz

#phased filtered 1000 genomes
kg=/local/chib/toconnor_grp/dloesch/1000_GENOMES/CCDG_14151_B01_GRM_WGS_2020-08-05_chr4.filtered.shapeit2-duohmm-phased.vcf.gz
bcftools view -r chr4:${start}-${end} -m2 -M2 -v snps -Oz -o 1kg.SNCA.vcf.gz --threads 4 $kg
#remove multiallelic sites and filter
bcftools view -e ID=@./data/1kg.SNCA.multiallelic_drop_freq.txt -g^miss -m2 -M2 -v snps -c 1 -o 1kg.SNCA.filtered.vcf.gz -Oz --threads 4 data/1kg.SNCA.vcf.gz 


#LARGE-PD
large=/local/chib/toconnor_grp/LARGE-PD/imputed_2021/chr4.dose.vcf.gz
bcftools view -r chr4:${start}-${end} -m2 -M2 -v snps -Oz -o large.SNCA.vcf.gz --threads 4 $large

#PGP
pgp=/local/chib/toconnor_grp/dloesch/peru/peru.chr4.HG38.vcf.gz
bcftools view -r chr4:${start}-${end} -m2 -M2 -v snps -Oz -o pgp.SNCA.vcf.gz --threads 4 $pgp

#TB
tb=/local/chib/oconnor_genomes/GLAD/data_after_imputation/luoetal_TBprogression/chr4.dose.vcf.gz
bcftools view -r chr4:${start}-${end} -m2 -M2 -v snps -Oz -o tb.SNCA.vcf.gz --threads 4 $tb

#subset regions

Rscript scripts/subset.R -s nalls.PD.sumstats.txt --geno large-PD.files.txt --prefix large.region --build hg38 -r 250
Rscript scripts/subset.R -s nalls.PD.sumstats.txt --geno 1KG.files.txt --prefix 1kg.region --build hg38 -r 250
