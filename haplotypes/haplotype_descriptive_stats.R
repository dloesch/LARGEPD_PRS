##script for getting descriptive statistics about SNCA haplotypes

popdata <- read.table("PD.merged_pops.pgp.SNCA.txt", header=FALSE, stringsAsFactors = FALSE)

pops <- c(unique(popdata$V3), unique(popdata$V2))
pops <- pops[!is.na(pops)]
pos <- 89704960

results <- as.data.frame(pops)

for(p in pops){
  
  block.file <- paste0("results/hap_blocks/SNCA/", p,".blocks.det")
  if(file.exists(block.file)){
    data <- data.table::fread(block.file, data.table=FALSE)
    data <- data[data$BP1 <= pos & data$BP2 >= pos,]
    results$KB[results$pops == p] <- ifelse(nrow(data) == 0, 0, data$KB)
    results$NSNP[results$pops == p] <- ifelse(nrow(data) == 0, 0, data$NSNPS)
    results$rs356182[results$pops == p] <- ifelse(nrow(data) == 0, NA, 
                                                  ifelse(grep(pos, data$SNPS) == 1, TRUE, FALSE))
  }else{
    next
  }

}

write.table(results, "./results/hap_blocks/SNCA/SNCA.summary.txt", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)

kg <- read.table("integrated_call_samples_v3.20130502.ALL.panel", header=TRUE, stringsAsFactors = FALSE)
snps <- blocks$SNP
pops <- unique(kg$super_pop)
results <- as.data.frame(snps)

for(i in 1:nrow(blocks)){
  for(p in pops){
    pop <- unique(kg$pop[kg$super_pop == p])
    pop <- paste0("sizekb_",pop)
    foo <- blocks[c("SNP", pop)]
    foo <- unlist(unname(foo[i,2:ncol(foo)]))
    results[[paste0(p, "_MEAN")]][i] <- mean(foo)
    results[[paste0(p, "_SD")]][i] <- sd(foo)
    results[[paste0(p, "_RANGE")]][i] <- paste0(signif(min(foo), digits=4), "-", signif(max(foo), digits=4))
  }
}

for(p in pops){
  results[[p]] <- paste0(signif(results[[paste0(p, "_MEAN")]], digits = 4), " (", signif(results[[paste0(p, "_SD")]], digits=4), ")")
}

colnames(results)[1] <- "SNP"
results <- results[c("SNP", "AFR", "AFR_RANGE", "AMR", "AMR_RANGE", "EAS", "EAS_RANGE", "EUR", "EUR_RANGE","SAS", "SAS_RANGE")]

write.table(results, "results/hap_blocks/PDhap_blocks.variation_summary.txt", sep='\t', quote=FALSE, row.names = FALSE, col.names = TRUE)

