#!/usr/bin/Rscript

library(argparser)
library(data.table)
library(tools)
source("./scripts/PRS_functions.R")

#set argument parser
argp <- arg_parser("claculate PRS using external GWAS data")

#reguired arguments
argp <- add_argument(argp, "--sumstats", help="GWAS summary statistics", type = "character")
argp <- add_argument(argp, "--geno", help="target genotype file or textfile with filenames if split by chrom", type = "character")

#optional arguments
argp <- add_argument(argp, "--prefix", help="output prefix", type = "character", default="output")
argp <- add_argument(argp, "--build", help="genome build of data", type = "character", default = "hg19")
argp <- add_argument(argp, "--CGAT", help="remove ambiguous sites", type = "logical", default = TRUE)
argp <- add_argument(argp, "--imputed", help="imputed data or genotype data", default=FALSE)

#optional file name columns
argp <- add_argument(argp, "--snp", help="column name for snp IDs", type = "character", default = NA)
argp <- add_argument(argp, "--beta", help="column name for betas",type = "character", default = NA)
argp <- add_argument(argp, "--pval", help="column name for p-values",type = "character", default = NA)
argp <- add_argument(argp, "--chrom", help="column name for bp positions", type = "character", default = NA)
argp <- add_argument(argp, "--pos", help="column name for bp positions", type = "character", default = NA)
argp <- add_argument(argp, "--ref", help="reference/other allele", type = "character", default = NA)
argp <- add_argument(argp, "--alt", help="alternate/effect allele", type = "character", default = NA)

#parse arguments
argsv <- parse_args(argp)


#summary stats and genotype data
geno.file <- argsv$geno
sumstats.file <- argsv$sumstats
build <- toupper(argsv$build)

#check input and exit
if(is.na(geno.file) | is.na(sumstats.file)){
  print(argp)
  stop("Error: Please provide GWAS summary stats and target genotype file(s)")
}

data <- fread(sumstats.file, header=TRUE, data.table = FALSE)

#column names
snpid <- argsv$snp
if(is.na(snpid)){
  snpid <- colnames(data)[colnames(data) %in% c("SNP", "snp", "Snp", "variant", "VARIANT", "ID")]
  if(length(snpid) != 1){
    stop("unable to infer SNP ID column name; please provide")
  }
}
beta <- argsv$beta
if(is.na(beta)){
  beta <- colnames(data)[colnames(data) %in% c("beta", "BETA", "weight")]
  OR <- colnames(data)[colnames(data) %in% c( "odds_ratio", "OR")]
  if(length(OR) == 1){
    beta <- "BETA"
    data$BETA <- log(data[[OR]])
  }
  if(length(beta) != 1){
    stop("unable to infer beta column name; please provide")
  }
}

pval <- argsv$pval
if(is.na(pval)){
  pval <- colnames(data)[colnames(data) %in% c("P", "PVALUE", "pvalue", "p", "Score.pval", "PVAL", "pval", "P_VALUE")]
  if(length(pval) != 1){
    stop("unable to infer pval column name; please provide")
  }
}

chrom <- argsv$chrom
if(is.na(chrom)){
  chrom <- colnames(data)[colnames(data) %in% c("C", "CHROM", "CHR", "c", "chrom", "chr", "chromosome", "CHROMOSOME")]
  if(length(chrom) != 1){
    stop("unable to infer chrom column name; please provide")
  }
}

pos <- argsv$pos
if(is.na(pos)){
  pos <- colnames(data)[colnames(data) %in% c("pos", "POS", "BP", "bp", "position", "POSITION")]
  if(length(pos) != 1){
    stop("unable to infer BP column name; please provide")
  }
}

ref <- argsv$ref
if(is.na(ref)){
  ref <- colnames(data)[colnames(data) %in% c("REF", "ref", "A2", "REFRENCE", "reference", "other_allele", "OTHER_ALLELE")]
  if(length(pos) != 1){
    stop("unable to infer ref/other allele column name; please provide")
  }
}

alt <- argsv$alt
if(is.na(alt)){
  alt <- colnames(data)[colnames(data) %in% c("alternate", "alt", "A1", "ALTERNATE", "alt", "ALT", "effect_allele", "EFFECT_ALLELE")]
  if(length(pos) != 1){
    stop("unable to infer alt/effect allele column name; please provide")
  }
}
#####parse sumstats file#######
#temporary prefix
temp <- paste0("temp", sample(1:1E6,1))

if(build == "HG38"){
  if(length(grep("chr",data[[chrom]])) == 0){
    data[[chrom]] <- paste0("chr", data[[chrom]])
  }
}


#prepare region file for subsetting genotype file
region <- data[c(chrom, pos)]
region.file <- paste0(temp,".region.txt")
write.table(region, region.file, sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)

#organize data
data <- data[c(snpid, chrom, pos, ref, alt, beta, pval)]
colnames(data) <- c("SNP", "CHROM", "POS", "REF", "ALT", "BETA", "PVAL")


####################################################
#infer if given single genotype file or file list
isVCF <- file_ext(sub(".gz", "", geno.file)) == "vcf"
isBCF <- file_ext(geno.file) == "bcf"
isBED <- file_ext(geno.file) == "bed"
isTXT <- file_ext(geno.file) == "txt"


##subset and read in genotype data
if(isVCF | isBCF){
  if(!file.exists(paste0(geno.file,".tbi"))){
    system(paste("tabix -f", geno.file))
  }
  outfile <- paste0(temp, ".vcf")
  system2("bcftools", 
          paste("view -R",region.file, "-Ov -o",outfile,geno.file))
  
  geno.data <- read.table(outfile, header=FALSE, colClasses = c("character"), comment.char="#")
}else if(isBED){
  #convert to vcf, if needed
  pfile <- sub(".bed", "", geno.file)
  tempvcf <- paste0(temp, ".2")
  system2("plink", 
          paste("--bfile", pfile, "--recode vcf-iid bgz  --out", tempvcf),
          stdout = FALSE)
  geno.file <- paste0(tempvcf, ".vcf.gz")
  #index
  if(!file.exists(paste0(geno.file,".tbi"))){
    system(paste("tabix -f", geno.file))
  }
  outfile <- paste0(temp, ".vcf")
  system2("bcftools", 
          paste("view -R",region.file, "-Ov -o",outfile,geno.file))
  
  geno.data <- read.table(outfile, header=FALSE, colClasses = c("character"), comment.char="#")
  
} else if(isTXT){
  ####folowing procedure is for a file list of genotype files####
  genolist <- fread(geno.file, header=FALSE, data.table=FALSE)
  if(ncol(genolist) < 1){
    stop("error: genotype filename list needs to be plain text, with one column and one filename per row")
  }
  
  if(nrow(genolist) < 22 | nrow(genolist) > 22){
    stop("error: do not see 22 files (one per autusomal chromosome)")
  }
  
  for(chr in 1:22){
    gfile <- genolist$V1[chr]
    outfile <- paste0(temp, ".chr", chr, ".vcf")
    
    #convert from plink to vcf, if needed
    isBED <- file_ext(gfile) == "bed"
    if(isBED){
      pfile <- sub(".bed", "", gfile)
      tempvcf <- paste0(temp, ".2.chr", chr)
      system2("plink", 
              paste("--bfile", pfile, "--recode vcf-iid bgz  --out", tempvcf),
              stdout = FALSE)
      gfile <- tempvcf
    }
    
    if(!file.exists(paste0(gfile,".tbi"))){
      system(paste("tabix -f", gfile))
    }
    
    #subset
    system2("bcftools", 
            paste("view -R",region.file, "-Ov -o",outfile,gfile))
    
    cat(paste0(outfile, "\n"), file=paste0(temp,".files.txt"),append = TRUE)

  }
  
  #concatenate
  outfile <- paste0(temp, ".merged.vcf")
  system2("bcftools", 
          paste("concat --file-list", paste0(temp,".files.txt"), "-Ov -o", outfile),
          stdout=FALSE)
  
  geno.data <- read.table(outfile, header=FALSE, colClasses = c("character"))
}else{
  stop("error: genotype file needs to be vcf, bcf, bed, or plain text file listing genotype files")
}


####process genotype data#####
geno.info <- geno.data[c(1:5)]
colnames(geno.info) <- c("CHROM", "POS", "SNP", "REF", "ALT")

#rename SNPs
geno.info$old_ids <- geno.info$SNP
geno.info$SNP <- paste(geno.info$CHROM, geno.info$POS, geno.info$REF, geno.info$ALT, sep=":")
geno.data[,3] <- geno.info$SNP

##check against sumstats file
foo <- sapply(geno.info$SNP, check_sumstats, simplify = FALSE)
info <- c()
for(i in 1:length(foo)){info <- c(info,foo[[i]][1])}
snps <- c()
for(i in 1:length(foo)){snps <- c(snps,foo[[i]][2])}

geno.info$INFO <- info
geno.info$SNP2 <- snps

#check if CG/AT
if(argsv$CGAT == TRUE){
  geno.info$CGAT <- mapply(find_CGAT, geno.info$REF, geno.info$ALT)
}else{
  geno.info$CGAT <- FALSE
}


#run filter
geno.info$FILTER <- mapply(set_filter,geno.info$INFO, geno.info$CGAT)
geno.info <- geno.info[geno.info$FILTER == "PASS",]

#write-out snps used in model
model_snps <- geno.info[c("old_ids", "SNP", "SNP2")]
colnames(model_snps) <- c("ORIG.ID", "SUMSTATS.ID")
write.table(model_snps, paste0(argsv$prefix, ".model_snps.txt"), sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

#filter geno data
geno.data <- geno.data[geno.data[,3] %in% geno.info$SNP,]
geno.data[,3] <- geno.info$SNP2

#set info
geno.data[,8] <- geno.info$INFO

#get subject ids from vcf
ids <- get_ids(input.vcf = outfile)

#convert to dosages or extract dosages if imputed
geno.dosages <- as.data.frame(ids)
gwas_weights <- data[c("SNP", "BETA")]

for(snp in geno.info$SNP2){
  if(argsv$imputed == TRUE){
    geno.dosages[[snp]] <- get_imputed_dosage(geno.data[geno.data[,3] == snp,])
  }else{
    #convert genotype calls to dosages
    geno.dosages[[snp]] <- to_dosage(geno.data[geno.data[,3] == snp,])
  }
  
  
  #flip direction of effect if necessary
  info <- geno.info$INFO[geno.info$SNP2 == snp]
  if(info == "ALLELE_SWITCH" | info == "ALLELE_SWITCH_STRAND_FLIP"){
    geno.dosages[[snp]] <- geno.dosages[[snp]]*data$BETA[data$SNP == snp]*-1
  }else{
    #multiply dosage by weight
    geno.dosages[[snp]] <- geno.dosages[[snp]]*data$BETA[data$SNP == snp]
  }

}

#sum weighted dosages to get PRS
geno.dosages$PRS <- rowSums(geno.dosages[,2:ncol(geno.dosages)], na.rm = TRUE)

#prepare PRS file
PRS <- geno.dosages[c("ids", "PRS")]
colnames(PRS) <- c("ID", "PRS_RAW")
PRS$PRS_AVG <- PRS$PRS_RAW/length(geno.info$SNP)
PRS$PRS_SCALED <- scale(PRS$PRS_RAW)

#name and write ouput
output.file <- paste0(argsv$prefix, ".PRS.txt")
write.table(PRS, output.file, col.names = TRUE, sep='\t', quote=FALSE, row.names = FALSE)

#cleanup
system2("mv", paste(outfile, paste0(argsv$prefix, ".PRS.vcf")))
system2("bgzip", paste0(argsv$prefix, ".PRS.vcf"))
system2("rm", paste0(temp, "*"))
