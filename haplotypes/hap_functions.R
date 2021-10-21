####subset vcf function ##

subset_vcf <- function(chr,pos,buffer){
  
  vcffile <- paste0("/mnt/hdd/1000_genomes/ALL.chr", 
                    chr,
                    ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz")
  
  start <- pos - buffer*1000 
  end <- pos + buffer*1000
  random <- sample(1:1000000000,1)
  tempvcf <- paste0("temp", random,".vcf.gz")
  
  #subset file
  system(paste("bcftools view -r",paste0(chr,":", start,"-", end)," -q 0.05 -m2 -M2 -v snps", vcffile, "| grep -v '##' >" , tempvcf))
  
  #read subsetted vcf 
  vcf <- read.table(tempvcf, header=TRUE, comment.char = "", stringsAsFactors = FALSE)
  
  #clean-up
  
  system(paste("rm ",tempvcf))
  
  return(vcf)
  
  
}

##subset vcf function version 2

subset_vcf2 <- function(vcffile,chr,start,end){
  
  
  random <- sample(1:1000000000,1)
  tempvcf <- paste0("temp", random,".vcf.gz")
  
  #subset file
  system(paste("bcftools view -r",paste0(chr,":", start,"-", end)," -q 0.05 -m2 -M2 -v snps", vcffile, "| grep -v '##' >" , tempvcf))
  
  #read subsetted vcf 
  vcf <- read.table(tempvcf, header=TRUE, comment.char = "", stringsAsFactors = FALSE)
  
  #clean-up
  
  system(paste("rm ",tempvcf))
  
  return(vcf)
  
  
}

#function for getting haplotypes from a phased vcf file, splitting genotypes at "|"
get_haps <- function(vcf){
  AA <- vcf[c(3:5)]
  
  snps <- vcf[,3]
  rownames(vcf) <- snps
  vcf <- vcf[c(10:ncol(vcf))]
  vcf <- as.data.frame(t(vcf))
  colnames(vcf) <- snps
  #define haplotypes
  hap1 <- as.data.frame(paste0(rownames(vcf), ".1"))
  colnames(hap1) <- "ID"
  hap2 <- as.data.frame(paste0(rownames(vcf), ".2"))
  colnames(hap2) <- "ID"
  
  snps <- colnames(vcf)
  
  get_A1 <- function(x){
    return(unlist(unname(strsplit(x, split="")))[1])
  }
  
  get_A2 <- function(x){
    return(unlist(unname(strsplit(x, split="")))[3])
  }
  
  for(snp in snps){
    vcf[[snp]] <- as.character(vcf[[snp]])
    vcf[[snp]] <- strtrim(vcf[[snp]],3)
    hap1[[snp]] <- sapply(X = vcf[[snp]], FUN = get_A1)
    hap1[[snp]] <- ifelse(hap1[[snp]] == 0, AA$REF[AA$ID == snp], AA$ALT[AA$ID == snp])
    
    hap2[[snp]] <- sapply(X = vcf[[snp]], FUN = get_A2)
    hap2[[snp]] <- ifelse(hap2[[snp]] == 0, AA$REF[AA$ID == snp], AA$ALT[AA$ID == snp])
  }
  
  
  for(i in 1:nrow(hap1)){
    
    hap1$haps[i] <- paste(unlist(unname(hap1[i,2:length(snps)])),collapse ="")
    
  }
  
  for(i in 1:nrow(hap2)){
    hap2$haps[i] <- paste(unlist(unname(hap2[i,2:length(snps)])),collapse ="")
  }
  
  haps <- rbind(hap1, hap2)
  
  return(haps)
}

#function for summarizing haplotypes when providing SNP (to keep track of REF/ALT status) and file with population data
summarize_haplodata <- function(vcf, haps, snp, pop.file){
  
  #read in population information
  pops <- read.table(pop.file, header=TRUE, stringsAsFactors = FALSE)
  colnames(pops) <- c("sample", "super_pop")
  #replace "-" with "." in IDS
  pops$sample <- gsub("-",".", pops$sample)
  pops$ID1 <- paste0(pops$sample, ".1")
  pops$ID2 <- paste0(pops$sample, ".2")
  
  #strip leading X from ID column
  haps$ID <- gsub("^X","",haps$ID)
  

  #keep intersection between popdata and haplodata
  haps$ID <- as.character(haps$ID)
  haps <- haps[haps$ID %in% c(pops$ID1, pops$ID2),]

  #format file
  haplotypes <- unique(haps$haps)
  
  haplodata <- as.data.frame(haplotypes)
  
  haplodata$hapID <- paste0("hap", 1:nrow(haplodata))
  
  #specify snp as target
  target <- snp
  AA <- vcf[c(3:5)]
  
  #format file
  haps <- merge(haps, haplodata, by.x="haps", by.y="haplotypes", all.x=TRUE)
  haps$target <- ifelse(haps[[target]] == AA$REF[AA$ID == target], "REF", "ALT")
  
  
  for(i in 1:nrow(haplodata)){
    haplodata$freq[i] <- nrow(haps[haps$hapID == haplodata$hapID[i],])/nrow(haps)
    
    for(p in unique(pops$super_pop)){
      samples <- c(pops$ID1[pops$super_pop == p], pops$ID2[pops$super_pop == p])
      haplodata[[p]][i] <- nrow(haps[haps$hapID == haplodata$hapID[i] & haps$ID %in% samples,])/nrow(haps[haps$ID %in% samples,])
    }
    
  }
  
  #get ref/alt status
  temp <- haps[c("hapID", "target")]
  temp <- temp[!duplicated(temp),]
  haplodata <- merge(haplodata, temp, by.x = "hapID", by.y="hapID", all.x=TRUE)
  
  #summarize data by pop
  
  data <- as.data.frame(unique(pops$super_pop))
  colnames(data) <- "pops"
  for(pop in unique(pops$super_pop)){
    samples <- c(pops$ID1[pops$super_pop == pop], pops$ID2[pops$super_pop == pop])
    N <- nrow(haps[haps$ID %in% samples,])
    ref <- nrow(haps[haps$ID %in% samples & haps$target == "REF",])
    alt <- nrow(haps[haps$ID %in% samples & haps$target == "ALT",])
    data$freqALT[data$pops == pop]  <- alt/N
    #ref
    topREF <- haplodata$hapID[haplodata[[pop]] == max(haplodata[[pop]][haplodata$target == "REF"])][1]
    data$topREF[data$pops == pop]  <- topREF
    data$topREF_freq[data$pops == pop] <- haplodata[[pop]][haplodata$hapID == topREF]
    data$propREF[data$pops == pop] <- nrow(haps[haps$hapID == topREF & haps$ID %in% samples,])/ref
    #alt
    topALT <- haplodata$hapID[haplodata[[pop]] == max(haplodata[[pop]][haplodata$target == "ALT"])][1]
    data$topALT[data$pops == pop]  <- topALT
    data$topALT_freq[data$pops == pop] <- haplodata[[pop]][haplodata$hapID == topALT]
    data$propALT[data$pops == pop] <- nrow(haps[haps$hapID == topALT & haps$ID %in% samples,])/alt
  }
  
  
  
  data$otherREF <- (data$topREF_freq/data$propREF) - data$topREF_freq
  data$otherALT <- (data$topALT_freq/data$propALT) - data$topALT_freq
  
  labels <- data[c("pops", "topREF", "topALT")]
  data <- data[c("pops", "topREF_freq", "otherREF", "topALT_freq", "otherALT")]
  
  haplodata <- list(haplodata, data, labels)
  
  return(haplodata)
}




#function for getting LD structure at a snp from 1000 Genomes 

get_LD <- function(snp, chr, pos){
  popfile <- "integrated_call_samples_v3.20130502.ALL.panel"
  
  start <- pos - 500000
  end <- pos + 500000
  
  random <- sample(1:1000000000,1)
  tempvcf <- paste0("temp", random,".vcf")
  
  vcffile <- paste0("/mnt/hdd/1000_genomes/ALL.chr", 
                    chr,
                    ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz")
  
  system(paste("bcftools view -r",paste0(chr,":", start,"-", end)," -q 0.05 -m2 -M2 -v snps", vcffile, "| grep -v '##' >" , tempvcf))
  
  for(pop in c("AFR", "AMR", "EAS", "EUR", "SAS")){
    system(paste("grep", pop, popfile, "| cut -f1 >", paste0(pop,".",random,".txt")))
    system(paste("plink --vcf", tempvcf, " --r2 inter-chr  --ld-snp", snp, "--out", paste0(pop,".",random), "--keep-fam", paste0(pop,".",random,".txt")))
  }
  
  #EUR
  file <- paste0("EUR.",random, ".ld" )
  EUR <- read.csv(file, sep="", stringsAsFactors=FALSE)
  EUR <- EUR[c(5:7)]
  colnames(EUR) <- c("EUR_BP","SNP", "EUR_R2")
  
  #AMR
  file <- paste0("AMR.",random, ".ld" )
  AMR <- read.csv(file, sep="", stringsAsFactors=FALSE)
  AMR <- AMR[c(5:7)]
  colnames(AMR) <- c("AMR_BP","SNP", "AMR_R2")
  
  data <- merge(EUR, AMR, by.x="SNP", by.y="SNP", all.x=TRUE, all.y="TRUE")
  
  
  #AFR
  file <- paste0("AFR.",random, ".ld" )
  AFR <- read.csv(file, sep="", stringsAsFactors=FALSE)
  AFR <- AFR[c(5:7)]
  colnames(AFR) <- c("AFR_BP","SNP", "AFR_R2")
  data <- merge(data, AFR, by.x="SNP", by.y="SNP", all.x=TRUE, all.y="TRUE")
  
  
  #EAS
  file <- paste0("EAS.",random, ".ld" )
  EAS <- read.csv(file, sep="", stringsAsFactors=FALSE)
  EAS <- EAS[c(5:7)]
  colnames(EAS) <- c("EAS_BP","SNP", "EAS_R2")
  data <- merge(data, EAS, by.x="SNP", by.y="SNP", all.x=TRUE, all.y="TRUE")
  
  
  #SAS
  file <- paste0("SAS.",random, ".ld" )
  SAS <- read.csv(file, sep="", stringsAsFactors=FALSE)
  SAS <- SAS[c(5:7)]
  colnames(SAS) <- c("SAS_BP","SNP", "SAS_R2")
  
  data <- merge(data, SAS, by.x="SNP", by.y="SNP", all.x=TRUE, all.y="TRUE")
  
  pops <- c("AFR", "AMR", "EAS", "EUR", "SAS")
  data$BP <- NA
  
  for(i in 1:nrow(data)){
    for(pop in pops){
      if(is.na(data$BP[i]) & !is.na(data[[paste0(pop,"_BP")]][i])){
        data$BP[i] <- data[[paste0(pop,"_BP")]][i]
      }else{
        next
      }
    }
  }
  
  data <- data[order(data$BP),]
  
  #cleanup
  files <- c(tempvcf)
  
  for(p in c("AFR", "AMR", "EAS", "EUR", "SAS")){
    files <- c(files, paste0(p, ".", random, ".*"))
  }
  files <- paste(files, collapse = " ")
  
  system(paste("rm", files))
  system(paste("rm", tempvcf))
  
  #return dataframe
  return(data)
}

#####get matrix of allele sharing ######

get_haplomat <- function(haplodata){
  top <- unique(c(haplodata[[3]]$topREF, haplodata[[3]]$topALT))
  top_haps <- haplodata[[1]][haplodata[[1]]$hapID %in% top,]
  top_haps$haplotypes <- as.character(top_haps$haplotypes)
  temp <- as.data.frame(1:nchar(top_haps$haplotypes[1])) 
  
  for(i in 1:nrow(top_haps)){
    temp[[top_haps$hapID[i]]] <- unlist(unname(strsplit(top_haps$haplotypes[i], split="")))
  }
  
  haplomat <- as.data.frame(matrix(ncol=length(top), nrow = length(top)))
  
  colnames(haplomat) <- top
  rownames(haplomat) <- top
  haplomat$ID <- top
  
  for(hap in top){
    for(h in top){
      h1 <- temp[[hap]]
      h2 <- temp[[h]]
      test <- h1 == h2
      N_same <- length(test[test == TRUE])/length(h1)
      haplomat[[h]][haplomat$ID == hap] <- N_same
    }
  }
  
  haplomat$ID <- NULL
  haplomat <- as.matrix(haplomat)
  
  return(haplomat)
}


####prep haplomat for plotting#######
prep_haplomat <- function(haplomat){
  require(reshape2)
  reorder_cormat <- function(cormat){
    # Use correlation between variables as distance
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat <-cormat[hc$order, hc$order]
  }
  
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  
  cormat <- reorder_cormat(haplomat)
  
  upper <- get_upper_tri(cormat)
  
  for(i in 1:ncol(upper)){
    upper[,i][upper[,i] == 1] <- NA
  }
  
  melted_cormat <- melt(upper, na.rm = TRUE)
  
  return(melted_cormat)
}
