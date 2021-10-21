
flip_allele <- function(x){
  y <- ifelse(x == "A", "T", ifelse(x == "T", "A", ifelse(x == "C", "G", ifelse(x == "G", "C", NA))))
  return(y)
}



##get dosages######
to_dosage <- function(gt){
  v <- unlist(unname(gt[1,10:ncol(gt)]))
  
  #if imputed, extracts genotype hard codes
  if(length(v[grep(":",v)]) > 0){v <- substr(v, 1, regexpr("\\:", v)-1)}
  
  v <- gsub("1.1",2,v)
  v <- gsub("0.1",1,v)
  v <- gsub("1.0",1,v)
  v <- gsub("0.0",0,v)
  v <- gsub("...", NA, v)
  v <- as.numeric(v)
  
  return(v)
}



####for imputed data####

get_imputed_dosage <- function(gt){
  v <- unlist(unname(gt[1,10:ncol(gt)]))
  foo <- c()
  for(i in 1:length(v)){
    foo <- c(foo, unlist(unname(strsplit(v[i], split=":")))[2])
  }
  foo <- as.numeric(foo)
  return(foo)
}



#################
check_sumstats <- function(target.snp){
  
  target.data <- geno.info[geno.info$SNP == target.snp,]
  chr <- target.data$CHROM
  pos <- target.data$POS
  ref <- target.data$REF
  alt <- target.data$ALT
  
  
  
  labels <- c("NOT_IN_REFERENCE", "SAME_STRAND", "ALLELE_SWITCH", "STRAND_FLIP", "STRAND_FLIP_ALLELE_SWITCH", "ALLELE_MISMATCH")
  PASS <- c( "SAME_STRAND", "ALLELE_SWITCH", "STRAND_FLIP", "STRAND_FLIP_ALLELE_SWITCH")
  
  sumstats <- data
  ref.data <- sumstats[sumstats$CHROM == chr & sumstats$POS == pos,]

  r.snp <- ref.data$SNP
  r.ref <- ref.data$REF
  r.alt <- ref.data$ALT
      
      
  label <- ifelse(r.ref == ref & r.alt == alt, labels[2], 
                      ifelse(r.ref == alt & r.alt == ref, labels[3],
                             ifelse(r.ref == flip_allele(ref) & r.alt == flip_allele(alt), labels[4],
                                    ifelse(r.ref == flip_allele(alt) & r.alt == flip_allele(ref), labels[5], labels[6]))))
      
  filter <- ifelse(label %in% PASS, "PASS", "FAIL")
  
  return(c(label, r.snp))
  
}


######
find_CGAT <- function(ref,alt){
  CGAT_filter <- ifelse(ref == "C" & alt == "G", TRUE, ifelse(ref == "G" & alt == "C", TRUE,
                                                                 ifelse(ref == "A" & alt == "T", TRUE,
                                                                        ifelse(ref == "T" & alt == "A", TRUE, FALSE))))
  return(CGAT_filter)
}
#####

set_filter <- function(info, CGAT){
  PASS <- c("SAME_STRAND", "STRAND_FLIP", "ALLELE_SWITCH", "STRAND_FLIP_ALLELE_SWITCH")
  filter <- ifelse(info %in% PASS, "PASS", "FAIL")
  
  filter <- ifelse(CGAT == TRUE, "FAIL", filter)
  
  return(filter)
}


###get subject ids from a vcf####


get_ids <- function(input.vcf, gz=FALSE){
  
  #get orginal header
  if(gz == TRUE){
    f <- gzfile(input.vcf, open="r")
  }else{
    f <- file(input.vcf, open="r")
  }
  
  header <- c()
  while(TRUE) {
    line = readLines(f, 1)
    if(length(line) == 0) break
    else if(grepl("^\\s*#{1}", line)) header <- c(header, line)
    else (break)
  }
  close(f)
  
  ids <- header[length(header)]
  ids <- unlist(unname(strsplit(ids, split="\t")))
  ids <- ids[10:length(ids)]
  
  return(ids)
}