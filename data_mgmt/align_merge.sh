#!/usr/bin/bash
#$ -q threaded.q  
#$ -P toconnor-lab
#$ -l mem_free=10G
#$ -e align.log
#$ -o align.log
#$ -pe thread 8
#$ -cwd
#$ -N align


###plink
export PATH=/usr/local/packages/plink-1.90.beta-3.6/bin:$PATH


###run align.R for merging ###
#align against 1KG
ref=/local/projects-t3/toconnor_grp/dloesch/PRS/LARGE-PD/data/1kg.SNCA.filtered.vcf.gz

#large-PD
target=/local/projects-t3/toconnor_grp/dloesch/PRS/LARGE-PD/data/large.SNCA.vcf.gz
prefix=large.SNCA

Rscript /local/projects-t3/toconnor_grp/dloesch/PRS/LARGE-PD/scripts/align/align.R -t $target -r $ref -p $prefix --threads 8 -C FALSE --build hg38

#pgp
target=/local/projects-t3/toconnor_grp/dloesch/PRS/LARGE-PD/data/pgp.SNCA.renamed.vcf.gz
prefix=pgp.SNCA      

Rscript /local/projects-t3/toconnor_grp/dloesch/PRS/LARGE-PD/scripts/align/align.R -t $target -r $ref -p $prefix --threads 8 -C FALSE --build hg38


#ipdgc
target=/local/projects-t3/toconnor_grp/dloesch/PRS/LARGE-PD/data/EUR.SNCA.vcf.gz
prefix=ipdgc.SNCA      

Rscript /local/projects-t3/toconnor_grp/dloesch/PRS/LARGE-PD/scripts/align/align.R -t $target -r $ref -p $prefix --threads 8 -C FALSE --build hg38

#tb cohort
target=/local/projects-t3/toconnor_grp/dloesch/PRS/LARGE-PD/data/tb.SNCA.vcf.gz
prefix=tb.SNCA

Rscript /local/projects-t3/toconnor_grp/dloesch/PRS/LARGE-PD/scripts/align/align.R -t $target -r $ref -p $prefix --threads 8 -C FALSE --build hg38

###merge####
#rename samples in ipdgc data, plink couldn't handle underscores
bcftools reheader -o ipdgc.SNCA.align.vcf.gz -s ids.txt ./data/ipdgc.SNCA.align.vcf.gz --threads 2 
mv ipdgc.SNCA.align.vcf.gz ./data #replace 
tabix -f ./data/ipdgc.SNCA.align.vcf.gz #re-index

#merge aligned files
bcftools merge -m snps -Oz -o temp.merged.vcf.gz ./data/ipdgc.SNCA.align.vcf.gz ./data/large.SNCA.align.vcf.gz ./data/tb.SNCA.align.vcf.gz ./data/pgp.SNCA.align.vcf.gz --threads 2
bcftools view -i 'F_MISSING<0.1' -Oz -o ipdgc_large_tb_pgp.merged.vcf.gz temp.merged.vcf.gz --threads 2

