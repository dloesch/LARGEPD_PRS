#!/bin/bash
#$ -N phase
#$ -q all.q  
#$ -P toconnor-lab
#$ -l mem_free=50G
#$ -e PHASE.logs
#$ -o PHASE.logs
#$ -cwd

chr=4
gt=$1
prefix=$2

beagle=/usr/local/packages/beagle-5.0/beagle.jar

out=$prefix.merged.chr$chr.PHASED
map=/local/chib/toconnor_grp/LARGE-PD/genetic_maps/plink.chr$chr.v2.GRCh38.map

ref=/local/projects-t3/toconnor_grp/dloesch/PRS/LARGE-PD/data/1kg.SNCA.filtered.vcf.gz

java -Xmx50g -jar $beagle gt=$gt out=$out map=$map chrom=chr$chr impute=false ref=$ref
