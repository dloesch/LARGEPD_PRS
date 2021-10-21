#!/usr/bin/Rscript

#path on cluster /usr/local/bin/Rscript
#######FOR ALIGNING GENOTYPE DATA AGAINST IMPUTATION REFERENCE ################

start_time <- Sys.time()

library(data.table)
library(parallel)
library(argparser)
library(tools)

#load supporting functions
source("/local/projects-t3/toconnor_grp/dloesch/PRS/LARGE-PD/scripts/align/align_functions.R")
source("/local/projects-t3/toconnor_grp/dloesch/PRS/LARGE-PD/scripts/align/align_utils.R")

#set argument parser
argp <- arg_parser("Align VCF file against a reference")

#reguired arguments
argp <- add_argument(argp, "--target", help="target vcf file", type = "character")
argp <- add_argument(argp, "--ref", help="reference bcf/vcf file", type = "character")

#optional arguments
argp <- add_argument(argp, "--prefix", help="output prefix", type = "character", default="output")
argp <- add_argument(argp, "--threads", help="number of threads", default = 2)
argp <- add_argument(argp, "--CGAT", help="TRUE: tests CG/AT sites, FALSE: removes CG/AT sites", type = "logical", default = FALSE)
argp <- add_argument(argp, "--filter", help="TRUE: filters output", type = "logical", default = FALSE)
argp <- add_argument(argp, "--merge", help="TRUE: merges output with reference", type = "logical", default = FALSE)
argp <- add_argument(argp, "--prune", help="TRUE: performs LD pruning", type = "logical", default = FALSE)
argp <- add_argument(argp, "--strict", help="TRUE: requires both LD and AF evidence for CGAT sites; FALSE: requires only AF", default=TRUE)
argp <- add_argument(argp, "--build", help="genome build: hg19 or hg38", type="character", default="hg19")

#parse arguments
argsv <- parse_args(argp)

#input parameters

input.vcf <- argsv$target
ref.file <- argsv$ref
out.prefix <- argsv$prefix
build <- toupper(argsv$build)


#additional parameters
test_CGAT <- as.logical(argsv$CGAT)
multithread <- ifelse(is.na(argsv$threads), FALSE, TRUE)
n.cores <- as.integer(argsv$threads)
plink <- as.logical(argsv$plink)
filter_vcf <- as.logical(argsv$filter)
merge_vcf <- as.logical(argsv$merge)
prune_vcf <- as.logical(argsv$prune)

#check input and exit
if(is.na(ref.file) | is.na(input.vcf)){
  print(argp)
  stop("Please provide reference and target files")
}


#make sure vcf file is indexed
#index
if(!file.exists(paste0(input.vcf,".tbi"))){
  system(paste("tabix -f",input.vcf))
}

#make sure ref file is indexed
if(!file.exists(paste0(ref.file,".tbi"))){
  system(paste("tabix -f",ref.file))
}


#set up cluster
if(multithread == TRUE){
  cl <- makeCluster(getOption("cl.cores", n.cores))
}


#get relevant data from vcf file
print("Subsetting files")


##using bcftools
system2("bcftools",paste("query -f '%CHROM  %POS %ID %REF  %ALT{0}\n'",input.vcf, ">",  paste0(out.prefix,".txt")))

data <- fread(file=paste0(out.prefix,".txt"), header=FALSE, stringsAsFactors = FALSE, colClasses = c("character"), data.table = FALSE)
data[,2] <- as.numeric(data[,2])

#get relevant data from ref file matched by position
#fwrite(data[c(1,2)], file=paste0(out.prefix,".region.txt"), sep='\t', col.names = FALSE)

system2("bcftools",paste("query -f '%CHROM  %POS %ID %REF  %ALT{0}\n'", ref.file, ">",  paste0(out.prefix,".ref.txt")))

ref.data <- fread(file=paste0(out.prefix,".ref.txt"), header=FALSE, stringsAsFactors = FALSE, colClasses = c("character"), data.table= FALSE)
ref.data[,2] <- as.numeric(ref.data[,2])
ref.data <- ref.data[ref.data$V2 %in% data$V2,]

#clean-up
system2("rm", paste0(out.prefix,".txt"))
#system2("rm", paste0(out.prefix,".region.txt"))
system2("rm", paste0(out.prefix,".ref.txt"))

print("Checking against reference")

checkstart <- Sys.time()
#run check alleles
if(multithread == TRUE){
  clusterExport(cl=cl, "ref.data")
  temp <- clusterMap(cl, check_alleles,data$V1, data$V2, data$V3, data$V4, data$V5, SIMPLIFY = FALSE)
  df <- data.frame(matrix(unlist(temp), nrow=length(temp), byrow=TRUE), stringsAsFactors = FALSE)
  
  data <- cbind(data, df)

}else{
  #version2, without multhreading:
  temp <- t(mapply(check_alleles,data$V1, data$V2, data$V3, data$V4, data$V5))
  data <- cbind(data, temp)
  
}
checkend <- Sys.time()
print("Check against ref runtime:")
print(checkend - checkstart)

#check if CG/AT
data$CGAT <- mapply(find_CGAT, data$V4, data$V5)

#update column names of logfile

colnames(data) <- c("CHROM", "POS", "ID", "REF", "ALT", "INFO", "SNP", "CGAT")

#flag duplicates
data$INFO <- as.character(data$INFO)
data$INFO <- ifelse(!is.na(data$SNP) & duplicated(data$SNP), "DUPLICATE", data$INFO)


####test CG/AT sites#####
data$POS <- as.numeric(data$POS)
if(test_CGAT == TRUE){
  print("Checking CG/AT sites")
  CGAT_start <- Sys.time()
  data$CGAT_filter <- NA
  if(multithread == TRUE){
    clusterExport(cl, "data")
    test <- clusterMap(cl, check_CGAT, snp=data$ID[data$CGAT == TRUE], info=data$INFO[data$CGAT == TRUE], 
                       target.file=input.vcf, ref.file = ref.file, strict=argsv$strict, prefix=out.prefix)
    data$CGAT_filter[data$CGAT == TRUE] <- unlist(unname(test))
      
  }else{
    test <- mapply(check_CGAT, snp=data$ID[data$CGAT == TRUE], info=data$INFO[data$CGAT == TRUE], 
                   target.file=input.vcf, ref.file = ref.file, strict=argsv$strict, prefix=out.prefix)
    data$CGAT_filter[data$CGAT == TRUE] <- unlist(unname(test))
  }
  
  CGAT_end <- Sys.time()
  #print run-time
  print("CGAT runtime:")
  print(CGAT_end - CGAT_start)
}else{
  data$CGAT_filter <- NA
}

#####apply filters#####

data$FILTER <- mapply(set_filter,data$INFO, data$CGAT_filter)

#make sure CG/AT sites that pass are flipped
data$INFO <- ifelse(!is.na(data$CGAT_filter) & data$CGAT_filter == "POSSIBLE_STRAND_FLIP" & data$INFO == "SAME_STRAND", 
                    "STRAND_FLIP_ALLELE_SWITCH", data$INFO)

data$INFO <- ifelse(!is.na(data$CGAT_filter) & data$CGAT_filter == "POSSIBLE_STRAND_FLIP" & data$INFO == "ALLELE_SWITCH", 
                    "STRAND_FLIP", data$INFO)

if(test_CGAT == FALSE){
  data$FILTER <- ifelse(data$CGAT == TRUE, "FAIL", data$FILTER)
}

#prepare data for output
log <- data
data <- data[data$FILTER != "FAIL",]

####prepare output vcf#####

print("Using PLINK to write new VCF")


data$REF2 <- ifelse(data$INFO == "STRAND_FLIP", flip_allele(data$REF), data$REF)
data$ALT2 <- ifelse(data$INFO == "STRAND_FLIP", flip_allele(data$ALT), data$ALT)

data$REF2 <- ifelse(data$INFO == "ALLELE_SWITCH", data$ALT, data$REF2)
data$ALT2 <- ifelse(data$INFO == "ALLELE_SWITCH", data$REF, data$ALT2)


data$REF2 <- ifelse(data$INFO == "STRAND_FLIP_ALLELE_SWITCH", flip_allele(data$ALT), data$REF2)
data$ALT2 <- ifelse(data$INFO == "STRAND_FLIP_ALLELE_SWITCH", flip_allele(data$REF), data$ALT2)


#write out update name
drop.file <- paste0(out.prefix,".temp.drop.txt")
fwrite(as.data.frame(unique(log$ID[log$FILTER == "FAIL"])), file=drop.file, sep='\t', quote=FALSE, col.names = FALSE)

name.file <- paste0(out.prefix,".temp.names.txt")
fwrite(data[c("ID", "SNP")], file=name.file, sep='\t', quote=FALSE, col.names = FALSE)

flip.file <- paste0(out.prefix,".temp.flip.txt")
temp <- data[data$INFO == "STRAND_FLIP" | data$INFO == "STRAND_FLIP_ALLELE_SWITCH",]
fwrite(temp[c("ID")], file=flip.file, sep='\t', quote=FALSE, col.names = FALSE)

a1.file <- paste0(out.prefix,".temp.a1.txt")
fwrite(data[c("SNP", "ALT2")], file=a1.file, sep='\t', quote=FALSE, col.names = FALSE)

plink.file <- paste0(out.prefix,".temp")
isBCF <- file_ext(input.vcf) == "bcf"
#drop failed variants and convert to plink
if(isBCF){
  system2("plink", paste("--bcf", input.vcf,"--exclude",drop.file, "--make-bed --out", plink.file),
          stdout = FALSE)
}else{
  system2("plink", 
          paste("--vcf", input.vcf,"--exclude",drop.file, "--make-bed --out", plink.file),
          stdout = FALSE)
}

#flip variants
plink.file2 <-  paste0(out.prefix,".temp2")
system2("plink", paste("--bfile", plink.file,"--flip", flip.file,  "--make-bed --out", plink.file2),
        stdout = FALSE)

final.file <- paste0(out.prefix, ".align")
chr.label <- ifelse(build == "HG38", "--output-chr chr26", "--output-chr 26")
system2("plink", paste("--bfile", plink.file2, "--update-name", name.file,chr.label,"--a1-allele", a1.file, "2 1 --recode vcf-iid bgz --out", final.file),
        stdout = FALSE)

system2("rm", paste0(out.prefix, ".temp*"))


###write-out log file
logfile <- paste0(out.prefix, ".align.log")
fwrite(log, file=logfile, quote=FALSE,sep="\t",na = "NA")

end_time <- Sys.time()

print("overall run time:")
print(end_time - start_time)


######optional actions###########
if(filter_vcf == TRUE){
  target.file=paste0(out.prefix, ".vcf.gz")
  system2("plink", 
          paste("--vcf", target.file, chr.label,"--keep-allele-order --geno 0.1 --mind 0.1 --hwe 1E-06  --recode vcf-iid bgz --out", paste0(out.prefix, ".filtered")),
          stdout = FALSE)
  system2("rm", "*.nosex")
  system2("rm", paste0(output.prefix,".drop.txt"))
}

if(merge_vcf == TRUE){
  target.file=paste0(out.prefix, ".vcf.gz")
  n.cores <- ifelse(is.na(n.cores), 1, n.cores)
  system2("tabix", paste0(target.file))
  system2("bcftools", paste("merge -m none --threads", n.cores, "-Oz -o", paste0(out.prefix, ".temp.vcf.gz"), target.file, ref.file))
  system2("bcftools", paste("view -g ^miss -m2 -M2 -v snps --threads", n.cores, "-Oz -o", paste0(out.prefix, ".merged.vcf.gz"), paste0(out.prefix, ".temp.vcf.gz")))
  system2("rm", paste0(out.prefix, ".temp.vcf.gz"))
  
}

if(prune_vcf == TRUE){
  target.file=paste0(out.prefix, ".vcf.gz")
  system2("plink", 
          paste("--vcf", target.file, "--indep-pairwise 50 5 0.2 --out", paste0(out.prefix)),
          stdout = FALSE)
  
  system2("plink", 
          paste("--vcf", target.file, "--keep-allele-order --exclude", paste0(out.prefix, ".prune.out"),"--recode vcf-iid bgz --out", paste0(out.prefix, ".pruned")),
          stdout = FALSE)
  system2("rm", "*.nosex")
}
