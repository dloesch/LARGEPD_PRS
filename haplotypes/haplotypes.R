###Script for getting SNCA haplotypes

#load data.table and haplotype functions
library(data.table)
source("./scripts/hap_functions.R")

#prep population file#
#this file specifies population membership across all dtasets
#read in 1KG file
kg <- read.table("integrated_call_samples_v3.20130502.ALL.panel", header=TRUE, stringsAsFactors = FALSE)
kg <- kg[c("sample", "super_pop")]
colnames(kg) <- c("ID", "POP")

#read in large-PDt file and format data
large <- read.table("large.pheno.03_2021.txt", header=TRUE, stringsAsFactors = FALSE)
large$POP <- ifelse(large$PD_STATUS == 1, "LARGEPD_CASE", "LARGEPD_CONTROL")
large <- large[c("sample.id", "POP")]
colnames(large)[1] <- "ID"

#read in large-PD+TB joint file and format data
large_tb <- read.table("large_tb.merged.pheno.txt", header=TRUE, stringsAsFactors = FALSE)
large_tb <- large_tb[c("sample.id", "PD_STATUS", "STUDY")]
large_tb$POP <- ifelse(large_tb$PD_STATUS == 1, "LARGEPD_CASE", "LARGEPD_CONTROL")
large_tb$POP <- ifelse(large_tb$STUDY == "TB", "PERU_TB", large_tb$POP)
colnames(large_tb)[1] <- "ID"
large_tb$STUDY <- NULL
large_tb$PD_STATUS <- NULL

#read in IPDGC file and format data
ipdgc <- read.table("./data/PHENO_NACHO.txt", header=TRUE, stringsAsFactors = FALSE)
ipdgc <- ipdgc[c("IID", "STATUS")]
colnames(ipdgc) <- c("ID", "POP")
ipdgc$POP <- ifelse(ipdgc$POP == 2, "IPDGC_CASE", "IPDGC_CONTROL")

#Read in PGP ids and format data
pgp <- fread("./data/pgp.SNCA.align.vcf.gz", header=TRUE, data.table=FALSE)
pgp <- colnames(pgp)[10:ncol(pgp)]
labels <- gsub("-...","",pgp)
amerindian <- c("MOCHES", "MATZES","CHOPCCAS", "UROS")

pgp <- as.data.frame(pgp)
colnames(pgp) <- "ID"
pgp$POP <- labels
pgp$POP <- ifelse(pgp$POP %in% amerindian, "PGP_AMER", "PGP_MESTIZO")

#combine data
#pop.data <- rbind(kg, large_tb, ipdgc, pgp)
pop.data <- rbind(kg, large, ipdgc, pgp)
pop.data <- pop.data[complete.cases(pop.data),]
pop.file <- "hap.populations.txt"
write.table(pop.data, pop.file, sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)


#define region based on haplotype blocks inferred using plink
hap_block <- read.table("./results/hap_blocks/Peru.blocks.det", header=TRUE, stringsAsFactors = FALSE)

#hg38 position of rs356182
pos <- 89704960
hap_block <- hap_block[hap_block$BP1 <= pos & hap_block$BP2 >=pos,]
start <- hap_block$BP1
end <- hap_block$BP2
snps <- unlist(unname(strsplit(hap_block$SNPS, split="\\|")))

#subset vcf file
vcf.file <- "./data/1kg_ipdgc_large_pgp.phased.vcf.gz"
vcf <- subset_vcf2(vcffile = vcf.file, "chr4", start, end)

#restrict to SNPS in block
vcf <- vcf[vcf$ID %in% snps,]

haps <- get_haps(vcf)
snp <- "4:89704960:G:A"
#get haplotypes
haplodata <- summarize_haplodata(vcf, haps, snp=snp, pop.file = pop.file)

#save
save(haps, file="./results/haplotypes/haplotypes.RData")
save(haplodata, file="./results/haplotypes/haplodata.RData")

