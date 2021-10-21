##script for testing halotpyes for association with PD in LARGE-PD and IPDGC

prefix <- "large"

load("results/haplotypes/haplodata.RData")
load("results/haplotypes/haplotypes.RData")


haps$ID <- as.character(haps$ID)

ID <- unique(gsub("\\..$", "", haps$ID))
data <- as.data.frame(ID)
data$ID1 <- paste0(data$ID,".1")
data$ID2 <- paste0(data$ID,".2")

#filter out rare haplotypes
haplotypes <- haplodata[[1]]$hapID[haplodata[[1]]$freq > 0.001]
for(h in haplotypes){
  haplotype <- as.character(haplodata[[1]]$haplotypes[haplodata[[1]]$hapID == h])
  carriers <- haps$ID[haps$haps == haplotype]
  
  data[[h]] <- ifelse(data$ID1 %in% carriers & data$ID2 %in% carriers, 2, 
                      ifelse(data$ID1 %in% carriers | data$ID2 %in% carriers, 1, 0))
  
}

#read in large-PD data
pheno <- read.table("large.pheno.03_2021.txt", header=TRUE, stringsAsFactors = FALSE)
pheno <- read.table("large.10PC.pheno.txt", header=TRUE, stringsAsFactors = FALSE)
rel <- read.table("large_tb.unrelated_toberemoved.txt", header=FALSE)
pheno <- pheno[!pheno$sample.id %in% rel$V1,]

#remove X that R appended to the IDS, just in case
data$ID <- gsub("^X", "", data$ID)

pheno <- merge(pheno, data, by.x="sample.id", by.y="ID", all.x = TRUE, sort = FALSE)
pheno <- pheno[complete.cases(pheno$PD_STATUS),]

#filter haplotypes for a min frequency of 0.01
large_haps <- c()
for(h in haplotypes){
  freq <- sum(pheno[[h]])/(nrow(pheno)*2)
  if(freq < 0.01){
    pheno[[h]] <- NULL
  }else{
    large_haps <- c(large_haps, h)
  }
}

results <- as.data.frame(large_haps)
colnames(results) <- "hapID"

#descriptive test and test for association via logistic regression
for(h in large_haps){
  #frequencies
  results$rs356182[results$hapID == h] <- ifelse(haplodata[[1]]$target[haplodata[[1]]$hapID ==h] == "REF", "G","A")
  results$FREQ[results$hapID == h] <- sum(pheno[[h]])/(nrow(pheno)*2)
  results$FREQ_CASES[results$hapID == h] <- sum(pheno[[h]][pheno$PD_STATUS ==1])/(nrow(pheno[pheno$PD_STATUS ==1,])*2)
  results$FREQ_CONTROLS[results$hapID == h] <- sum(pheno[[h]][pheno$PD_STATUS ==0])/(nrow(pheno[pheno$PD_STATUS ==0,])*2)
  
  #logistic regression
  f=as.formula(paste0("PD_STATUS~AGE+SEX+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+SITE+",h))
  fit <- glm(f, family="binomial", data=pheno)
  pvals <- summary(fit)$coefficients
  
  #save results
  results$BETA[results$hapID == h] <- pvals[rownames(pvals) == h][1]
  results$SE[results$hapID == h] <- pvals[rownames(pvals) == h][2]
  results$PVAL[results$hapID == h] <-  pvals[rownames(pvals) == h][4]
  
}

#concodane of direction of effect with rs356182 allele status
results$CONCORDANCE <- ifelse(results$BETA < 0 & (results$FREQ_CASES < results$FREQ_CONTROLS), TRUE, 
                              ifelse(results$BETA > 0 & (results$FREQ_CASES > results$FREQ_CONTROLS), TRUE, FALSE))

#write-table
out.file <- paste0(prefix, ".rs356182_haplotypes.ASSOC.10PC.txt")
write.table(results, out.file, sep='\t', quote=FALSE, col.names = TRUE, row.names = FALSE)


#repeat with IPDGC data
prefix <- "IPDGC"
pheno <- read.table("data/PHENO_NACHO.txt", header = TRUE, stringsAsFactors = FALSE)

pheno$STATUS <- as.factor(pheno$STATUS)

pheno <- merge(pheno, data, by.x="IID", by.y="ID", all.x = TRUE, sort = FALSE)
pheno <- pheno[complete.cases(pheno),]

#filter haplotypes
ipdgc_haps <- c()
for(h in haplotypes){
  freq <- sum(pheno[[h]])/(nrow(pheno)*2)
  if(freq < 0.01){
    pheno[[h]] <- NULL
  }else{
    ipdgc_haps <- c(ipdgc_haps, h)
  }
}

results <- as.data.frame(ipdgc_haps)
colnames(results) <- "hapID"

for(h in ipdgc_haps){
  results$rs356182[results$hapID == h] <- ifelse(haplodata[[1]]$target[haplodata[[1]]$hapID ==h] == "REF", "G","A")
  results$FREQ[results$hapID == h] <- sum(pheno[[h]])/(nrow(pheno)*2)
  results$FREQ_CASES[results$hapID == h] <- sum(pheno[[h]][pheno$STATUS ==2])/(nrow(pheno[pheno$STATUS ==2,])*2)
  results$FREQ_CONTROLS[results$hapID == h] <- sum(pheno[[h]][pheno$STATUS ==1])/(nrow(pheno[pheno$STATUS ==1,])*2)
  
  f=as.formula(paste0("STATUS~AGE+SEX+PC1+PC2+PC3+PC4+PC5+COHORT+",h))
  fit <- glm(f, family="binomial", data=pheno)
  pvals <- summary(fit)$coefficients
  
  results$BETA[results$hapID == h] <- pvals[rownames(pvals) == h][1]
  results$SE[results$hapID == h] <- pvals[rownames(pvals) == h][2]
  results$PVAL[results$hapID == h] <-  pvals[rownames(pvals) == h][4]
  
}

results$CONCORDANCE <- ifelse(results$BETA < 0 & (results$FREQ_CASES < results$FREQ_CONTROLS), TRUE, 
                              ifelse(results$BETA > 0 & (results$FREQ_CASES > results$FREQ_CONTROLS), TRUE, FALSE))

out.file <- paste0(prefix, ".rs356182_haplotypes.ASSOC.txt")
write.table(results, out.file, sep='\t', quote=FALSE, col.names = TRUE, row.names = FALSE)

