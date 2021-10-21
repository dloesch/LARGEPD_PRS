#!/usr/bin/bash
#$ -q threaded.q
#$ -P toconnor-lab
#$ -pe thread 4
#$ -l mem_free=10G
#$ -e concat.log
#$ -o concat.log
#$ -cwd
#$ -N concat
#$ -hold_jid merge


#concatenate merged files from large_PD and tb study

echo -n > tb_large_files.txt
for chr in {1..22}; do echo tb_large.chr$chr.pruned.vcf.gz >> tb_large_files.txt; done

bcftools concat -f tb_large_files.txt -Oz -o large_tb.pruned.all.vcf.gz --threads 4
