#!/usr/local/bin/R


############check alleles function #######################
check_alleles <- function(chr, pos, snp, ref, alt){
  
  source("/local/projects-t3/toconnor_grp/dloesch/PRS/LARGE-PD/scripts/align/align_utils.R")
  pos <- as.numeric(pos)
  
  #subset ref data that is made available in environment
  r <- ref.data[ref.data[,2] == pos,]
  
  labels <- c("NOT_IN_REFERENCE", "SAME_STRAND", "ALLELE_SWITCH", "STRAND_FLIP", "STRAND_FLIP_ALLELE_SWITCH", "ALLELE_MISMATCH")
  PASS <- c( "SAME_STRAND", "ALLELE_SWITCH", "STRAND_FLIP", "STRAND_FLIP_ALLELE_SWITCH")
  
  if(nrow(r) == 0){
    label <- labels[1]
    filter <- "FAIL"
    r.snp <- NA
  }else{
    
    
    for(j in 1:nrow(r)){
      r.chr <- r[j,1]
      r.pos <- as.numeric(r[j,2])
      r.snp <- r[j,3]
      r.ref <- r[j,4]
      r.alt <- r[j,5]
      
      
      label <- ifelse(r.ref == ref & r.alt == alt, labels[2], 
                      ifelse(r.ref == alt & r.alt == ref, labels[3],
                             ifelse(r.ref == flip_allele(ref) & r.alt == flip_allele(alt), labels[4],
                                    ifelse(r.ref == flip_allele(alt) & r.alt == flip_allele(ref), labels[5], labels[6]))))
      
      filter <- ifelse(label %in% PASS, "PASS", "FAIL")
      
      #if label is assigned, break loop even if there is another snp at position
      if(filter != "FAIL"){
        break
      }
    }
  }
  
  return(c(label, r.snp))
  
}

##########find CGAT##############

find_CGAT <- function(ref,alt){
  CGAT_filter <- ifelse(ref == "C" & alt == "G", TRUE, ifelse(ref == "G" & alt == "C", TRUE,
                                                                            ifelse(ref == "A" & alt == "T", TRUE,
                                                                                   ifelse(ref == "T" & alt == "A", TRUE, FALSE))))
  return(CGAT_filter)
}


##check CGAT sites######
####check CGAT against reference######
check_CGAT <- function(snp,info,target.file, ref.file, strict=FALSE, prefix=NA){
  source("/local/projects-t3/toconnor_grp/dloesch/PRS/LARGE-PD/scripts/align/CGAT_tests.R")
  source("/local/projects-t3/toconnor_grp/dloesch/PRS/LARGE-PD/scripts/align/align_utils.R")
  
  #try hard-coding
  df <- data
  
  
  if(info == "NOT_IN_REFERENCE" | info == "DUPLICATE" | info == "ALLELE_MISMATCH"){
    AF_test <- NA
    R_test <- NA
    CGAT_filter <- NA
  }else{
    #subset target and ref files
    chr <- df$CHROM[df$ID == snp]
    pos <- df$POS[df$ID == snp]
    r.snp <- df$SNP[df$ID == snp]
    
    target.gt <- subset_vcf(target.file, chr=chr, pos=pos, snp=snp)
    ref.gt <- subset_vcf(ref.file, chr=chr, pos=pos, snp=r.snp)
    
    #check for allele switch
    allele_switch <- ifelse(info == "ALLELE_SWITCH", TRUE, FALSE)
    
    #evaluate tests
    
    if(strict == TRUE){
      #check allele frequencies, keep variants with MAF > 0.4 for now
      AF_test <- allelefreq_test(target.gt = target.gt, ref.gt = ref.gt, maf.filter=FALSE, allele_switch)
      #also requires correlation data
      #match sign of dosage correlation
      R_test <- r_flip_test(snp = snp, target.gt = target.gt, target.file=target.file,  ref.gt = ref.gt, ref.file = ref.file ,df = df)
      
      CGAT_filter <- ifelse(AF_test == "SAME_STRAND" & R_test == "SAME_STRAND", "SAME_STRAND", 
                            ifelse(AF_test == "POSSIBLE_STRAND_FLIP" & R_test == "POSSIBLE_STRAND_FLIP", "POSSIBLE_STRAND_FLIP", "UNINFORMATIVE"))
    }else{
      #check allele frequencies, filter out variants with MAF > 0.4
      AF_test <- allelefreq_test(target.gt = target.gt, ref.gt = ref.gt, maf.filter=TRUE, allele_switch)
      #just uses allele frequencies
      
      R_test <- NA
      CGAT_filter <- ifelse(AF_test == "SAME_STRAND", "SAME_STRAND", 
                            ifelse(AF_test == "POSSIBLE_STRAND_FLIP", "POSSIBLE_STRAND_FLIP", "UNINFORMATIVE"))
    }
    
  }
  
  prefix <- ifelse(is.na(prefix), "output", prefix)
  write(paste(snp, AF_test, R_test, CGAT_filter, sep =' '), file=paste0(prefix,".CGAT.log"), append = TRUE)
  return(CGAT_filter)
}


#####set filter#######
set_filter <- function(info, CGAT_filter=NA){
  PASS <- c("SAME_STRAND", "STRAND_FLIP", "ALLELE_SWITCH", "STRAND_FLIP_ALLELE_SWITCH")
  filter <- ifelse(info %in% PASS, "PASS", "FAIL")
  
  filter <- ifelse(!is.na(CGAT_filter) & (CGAT_filter == "SAME_STRAND" | CGAT_filter == "POSSIBLE_STRAND_FLIP") & info %in% PASS,"PASS",
                   ifelse(!is.na(CGAT_filter) & CGAT_filter == "UNINFORMATIVE","FAIL", filter))
  
  return(filter)
}

####write vcf function ######
###these functions are no longer used in this impelentation

write_vcf_header <- function(input.vcf, output.vcf){
  
  #get orginal header
  f <- gzfile(input.vcf, open="r")
  header <- c()
  while(TRUE) {
    line = readLines(f, 1)
    if(length(line) == 0) break
    else if(grepl("^\\s*#{1}", line)) header <- c(header, line)
    else (break)
  }
  close(f)
  
  #prepare output vcf
  output.header <- c("##fileformat=VCFv4.2", 
                     paste0("##filedate=", format(Sys.time(), "%Y%m%d")), 
                     "##source=align.R",
                     "##FORMAT=<ID=GT,Number=1,Type=String,Description='Genotype'",
                     header[length(header)])
  fileConn <- file(output.vcf)
  writeLines(output.header, fileConn)
  close(fileConn)
  
  print("vcf file initialized")
}

####base align_vcf function ######
align_vcf <- function(snp, info){
  
  source("/local/projects-t3/toconnor_grp/dloesch/PRS/LARGE-PD/scripts/align/align_utils.R")
  
  #hardcode
  #v <- vcf[vcf[,3] == snp,]
  v <- other[other[,3] == snp,]
  
  
  if( info == "SAME_STRAND"){
    return(v)
  }else if( info == "ALLELE_SWITCH"){
    ref <- v[1,4]
    alt <- v[1,5]
    v[1,4] <- alt
    v[1,5] <- ref
    
    #change 0 1 assignments
    v <- switch_vcf(v)
    
    return(v)
  }else if( info == "STRAND_FLIP"){
    v[1,4] <- flip_allele(v[1,4])
    v[1,5] <- flip_allele(v[1,5])
    
    return(v)
  }else if( info == "STRAND_FLIP_ALLELE_SWITCH"){
    
    #switch and flip alleles
    ref <- v[1,4]
    alt <- v[1,5]
    
    v[1,4] <- flip_allele(alt)
    v[1,5] <- flip_allele(ref)
    
    #change 0 1 assignments:
    v <- switch_vcf(v)
    return(v)
  }else{
    stop("error: only provide variants that pass filters to this function")
  }
}

