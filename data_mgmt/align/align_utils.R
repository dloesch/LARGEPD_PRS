#!/usr/local/bin/R

###align utilities######
###flip allele####

flip_allele <- function(x){
  y <- ifelse(x == "A", "T", ifelse(x == "T", "A", ifelse(x == "C", "G", ifelse(x == "G", "C", NA))))
  return(y)
}


#switch 1/0 assignment in vcf records

switch_vcf <- function(vcf){
  n <- ncol(vcf)
  temp <- t(vcf)
  temp <- as.data.frame(temp)
  temp[,1] <- as.character(temp[,1])
  temp[,1] <- gsub("\\|", "/", temp[,1])
  temp <- ifelse(temp[,1] == "1/1", "0/0", ifelse(temp[,1] == "0/1", "1/0", ifelse(temp[,1] == "0/0", "1/1", temp[,1])))
  temp <- as.character(temp)
  foo <- as.data.frame(matrix(ncol = ncol(vcf), nrow=1))
  foo[1,] <- temp
  return(foo)
}


##get dosages######
to_dosage <- function(gt){
  v <- unlist(unname(gt[1,10:ncol(gt)]))
  
  #if imputed, gets hard coded genotype
  if(length(v[grep(":",v)]) > 0){v <- substr(v, 1, regexpr("\\:", v)-1)}
  
  v <- gsub("1.1",2,v)
  v <- gsub("0.1",1,v)
  v <- gsub("1.0",1,v)
  v <- gsub("0.0",0,v)
  v <- gsub("...", NA, v)
  v <- as.numeric(v)
  
  return(v)
}

##flip dosages###
flip_dosage <- function(dosage){
  foo <-gsub(2,999, dosage)
  foo <- gsub(0,2, foo)
  foo <- gsub(999,0, foo)
  foo <- as.numeric(foo)
  return(foo)
}

##subset and read in vcf

subset_vcf <- function(input.vcf, chr, pos, snp=NA){
  
  require(data.table)
  
  foo <- sample(1:1E8,1)
  out.prefix <- paste0("temp.",foo)
  out.vcf <- paste0(out.prefix, ".vcf")
  system2("bcftools",paste("view -o", out.vcf, "-Ov -r", paste0(chr,":", pos), input.vcf))
  v <- fread(out.vcf, header = TRUE, data.table = FALSE)
  
  if(!is.na(snp)){
    v <- v[v[,3] == snp,]
  }
  
  #clean up
  system2("rm", paste0(out.prefix,"*" ))
  
  return(v)
}


