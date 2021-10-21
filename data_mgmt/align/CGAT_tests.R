
###allele frequency test#######

allelefreq_test <- function(target.gt,ref.gt, maf.filter= TRUE, allele_switch=FALSE){
  
source("/local/projects-t3/toconnor_grp/dloesch/PRS/LARGE-PD/scripts/align/align_utils.R")
  
 #get dosages 
  
  target <- to_dosage(target.gt)
  
  ref <- to_dosage(ref.gt)
  
  #calculate allele frequencies
  n.alleles <- length(target[!is.na(target)])*2
  af <- sum(target, na.rm=TRUE)/n.alleles
  
  r.alleles <- length(ref[!is.na(ref)])*2
  r.af <- sum(ref, na.rm = TRUE)/r.alleles  
  
  #flip af if allele switch
  if(allele_switch == TRUE ){
    af <- 1-af
  }else{
    af <- af
  }
  
  #check if abs diff > if "flipping"
  
  diff1 <- abs(af - r.af)
  diff2 <- abs((1-af) - r.af)

  AF_test <- ifelse(diff1 < diff2, "SAME_STRAND", "POSSIBLE_STRAND_FLIP")

  #filter if MAF is > 0.4 in one of the 2 datasets
  maf_test <- ifelse((af > 0.4 & af < 0.6) | (r.af > 0.4 & r.af < 0.6), "FAIL", "PASS")
  
  
  if(maf.filter== TRUE){
    AF_test <- ifelse(maf_test == "FAIL", "UNINFORMATIVE", AF_test)
  }
  
  return(AF_test)
}




#####allele correlations test#######

r_flip_test <- function(snp, target.gt, target.file, ref.gt, ref.file, df){
  
  source("/local/projects-t3/toconnor_grp/dloesch/PRS/LARGE-PD/scripts/align/align_utils.R")
  
  #check for allele switch
  allele_switch <- ifelse(df$INFO[df$ID == snp] == "ALLELE_SWITCH", TRUE, FALSE)
  
  #get adjacent snp
  df$DIFF <- abs(df$POS - df$POS[df$ID == snp])
  df$DIFF <- ifelse(df$DIFF == 0, NA, df$DIFF)
  df <- df[complete.cases(df$DIFF),]
  df <- df[order(df$DIFF),]
  
  #pick closest snp that is in reference and not an allele mismatch
  for(i in 1:nrow(df)){
    if(df$INFO[i] == "NOT_IN_REFERENCE" | df$INFO[i] == "ALLELE_MISMATCH"){
      next
    }else{
      snp2 <- df$ID[i]
      break
    }
  }
  
  #get 2nd reference snp name
  r.snp <- df$SNP[df$ID == snp2]
  
  #get second snp from vcf files
  chr <- df$CHROM[df$ID == snp2]
  pos <- df$POS[df$ID == snp2]
  
  target.gt2 <- subset_vcf(target.file, chr=chr, pos=pos, snp=snp2)
  ref.gt2 <- subset_vcf(ref.file, chr=chr, pos=pos, snp=r.snp)
  
  target <- to_dosage(target.gt)
  target2 <- to_dosage(target.gt2)
  
  if(allele_switch == TRUE){
    flipped <- target
    target <- flip_dosage(target)
  }else{
    flipped <- flip_dosage(target)
  }
  
  if(df$INFO[df$ID == snp2] == "ALLELE_SWITCH" | df$INFO[df$ID == snp2] == "STRAND_FLIP_ALLELE_SWITCH"){
    target2 <- flip_dosage(target2)
  }
  
  ref <- to_dosage(ref.gt)
  ref2 <- to_dosage(ref.gt2)
  
  target.cor <- cor(target, target2, use="complete.obs")
  
  flip.cor <- cor(flipped, target2, use = "complete.obs")
  ref.cor <- cor(ref, ref2, use="complete.obs")
  
  if(is.na(target.cor) | is.na(flip.cor) | is.na(ref.cor)){
    R_test <- "UNINFORMATIVE"
  }else{
    r1 <- ifelse(target.cor > 0 &  ref.cor > 0, "SAME", ifelse( target.cor < 0 & ref.cor < 0, "SAME", "OPPOSITE"))
    r2 <- ifelse(flip.cor > 0 &  ref.cor > 0, "SAME", ifelse(flip.cor < 0 & ref.cor < 0, "SAME", "OPPOSITE"))
    
    R_test <- ifelse(r1 == "SAME" & r2 == "OPPOSITE", "SAME_STRAND", ifelse(r1 == "OPPOSITE" & r2 == "SAME", "POSSIBLE_STRAND_FLIP", "UNINFORMATIVE"))
  }
  
  
  return(R_test)
}
