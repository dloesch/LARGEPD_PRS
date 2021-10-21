#!/usr/bin/bash
#$ -q all.q  
#$ -P toconnor-lab
#$ -l mem_free=5G
#$ -e blocks.log
#$ -o blocks.log
#$ -cwd
#$ -N blocks

##Script for estimating haplotype blocks using PLINK
###plink
export PATH=/usr/local/packages/plink-1.90.beta-3.6/bin:$PATH


#create merged population file
#1kg
awk -v OFS='\t' 'NR>1 {print $1, $2, $3}' integrated_call_samples_v3.20130502.ALL.panel > 1kg.txt
#largePD
awk -v OFS='\t' 'NR>1 {print $1, $5, $2}' large.pheno.03_2021.txt > large.txt
sed -i 's/\t0/\tLARGEPD_CONTROL/g' large.txt
sed -i 's/\t1/\tLARGEPD_CASE/g' large.txt
#ipdgc
awk -v OFS='\t' 'NR>1 {print $2, $6, $5}' ./data/PHENO_NACHO.txt > ipdgc.txt
sed -i 's/\t1/\tIPDGC_CONTROL/g' ipdgc.txt
sed -i 's/\t2/\tIPDGC_CASE/g' ipdgc.txt

#pgp
plink --vcf ./data/pgp.PRS.vcf.gz --make-just-fam --out temp.pgp --memory 600M
cut -f1 -d" " temp.pgp.fam > pgp.ids.txt
cut -f1 -d"-" temp.pgp.fam > pgp.pop.txt
paste pgp.ids.txt pgp.pop.txt > temp.pgp.txt
awk -v OFS='\t' '{print $1, $2, "PGP"}' temp.pgp.txt > pgp.txt

#1kg and large for PD regions
echo -e "ID\tPOP\tSUPERPOP" > PD.merged_pops.region.txt
cat 1kg.txt large.txt >> PD.merged_pops.region.txt

#ipdgc large and 1KG for SNCA
echo -e "ID\tPOP\tSUPERPOP" > PD.merged_pops.SNCA.txt
cat 1kg.txt large.txt ipdgc.txt >> PD.merged_pops.SNCA.txt

#ipdgc large 1KG 
echo -e "ID\tPOP\tSUPERPOP" > PD.merged_pops.pgp.SNCA.txt
cat 1kg.txt large.txt ipdgc.txt pgp.txt >> PD.merged_pops.pgp.SNCA.txt

#clean up temp files
rm 1kg.txt large.txt ipdgc.txt temp.pgp* pgp.ids.txt pgp.pop.txt

##estimating haplotype block size###
mkdir -p results/hap_blocks/SNCA

#convert vcf to plink file
#vcf=./data/1kg_large.merged.region.vcf.gz
vcf=./data/1kg_ipdgc_large_tb_pgp.merged.chr4.PHASED.vcf.gz
#subset
bcftools view -r chr4:89454960-89954960 -Oz -o temp.vcf.gz $vcf
vcf=temp.vcf.gz
plink --vcf $vcf --make-bed --out temp --memory 5000M

#superpops=$(cut -f3 PD.merged_pops.region.txt | sort -u | grep -v POP | grep -v NA)
superpops=$(cut -f3 PD.merged_pops.pgp.SNCA.txt | sort -u | grep -v POP | grep -v NA)

for pop in $superpops; do
	#cat PD.merged_pops.region.txt | grep $pop | cut -f1 > temp.txt
	cat PD.merged_pops.pgp.SNCA.txt | grep $pop | cut -f1 > temp.txt
	plink --bfile temp --maf 0.05 --hwe 1E-06 --geno 0.1 --mind 0.1 --blocks no-pheno-req --blocks-max-kb 10000 --out ./results/hap_blocks/SNCA/${pop} --keep-fam temp.txt --snps-only --memory 5000M
done

pops=$(cut -f2 PD.merged_pops.pgp.SNCA.txt | sort -u | grep -v POP | grep -v NA)

for pop in $pops; do
	#cat PD.merged_pops.region.txt | grep $pop | cut -f1 > temp.txt	
	cat PD.merged_pops.pgp.SNCA.txt | grep $pop | cut -f1 > temp.txt
	plink --bfile temp --maf 0.05 --hwe 1E-06 --geno 0.1 --mind 0.1 --blocks no-pheno-req --blocks-max-kb 10000 --out ./results/hap_blocks/SNCA/${pop} --keep-fam temp.txt --snps-only --memory 5000M
done
