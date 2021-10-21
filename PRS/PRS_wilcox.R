##script for characterizing PRS by LARGE-PD recruitment sites

#read in phenotype data
pheno <- read.table("large_tb.merged.pheno.txt", header=TRUE, stringsAsFactors = FALSE)

#read in PRS results
PRS <- read.delim("~/WORKSPACE/large.PRS.txt", stringsAsFactors=FALSE)

#merge PRS and phenotype data
data <- merge(pheno, PRS, by.x="sample.id", by.y="ID", all.x=TRUE, sort=FALSE)

results <- as.data.frame(sites)

data <- data[complete.cases(data),]

#fore each site, calculate mean and SD of PRS, plus case-control ratio
for(s in sites){
  results$CC_RATIO[results$sites == s] <- nrow(data[data$COUNTRY == s & data$PD_STATUS == 1,])/nrow(data[data$COUNTRY == s,])
  results$MEAN[results$sites == s] <- mean(data$PRS_RAW[data$COUNTRY == s])
  results$SD[results$sites == s] <- sd(data$PRS_RAW[data$COUNTRY == s])
  
  #using Peru as reference, test PRS distribution with Wilcox
  if(s == "Peru"){
    results$PVAL[results$sites == s] <- NA
  }else{
    test <- wilcox.test(data$PRS_RAW[data$COUNTRY == s], data$PRS_RAW[data$COUNTRY == "Peru"])
    results$PVAL[results$sites == s] <- test$p.value
  }
  
}

write.table(results, "large.PRS_summary.txt", sep='\t', quote=FALSE, col.names = TRUE, row.names = FALSE)

#now repeat using external controls from TB cohort:
PRS <- read.delim("~/WORKSPACE/large_tb.merged.PRS.txt", stringsAsFactors=FALSE)

data <- merge(pheno, PRS, by.x="sample.id", by.y="ID", all.x=TRUE, sort=FALSE)

results <- as.data.frame(sites)

data <- data[complete.cases(data),]

for(s in sites){
  results$CC_RATIO[results$sites == s] <- nrow(data[data$COUNTRY == s & data$PD_STATUS == 1,])/nrow(data[data$COUNTRY == s,])
  results$MEAN[results$sites == s] <- mean(data$PRS_RAW[data$COUNTRY == s])
  results$SD[results$sites == s] <- sd(data$PRS_RAW[data$COUNTRY == s])
  
  if(s == "Peru"){
    results$PVAL[results$sites == s] <- NA
  }else{
    test <- wilcox.test(data$PRS_RAW[data$COUNTRY == s], data$PRS_RAW[data$COUNTRY == "Peru"])
    results$PVAL[results$sites == s] <- test$p.value
  }
  
}

write.table(results, "large_tb.PRS_summary.txt", sep='\t', quote=FALSE, col.names = TRUE, row.names = FALSE)

#plot case-control ratio
library(ggplot2)

pdf("MEAN_PRS_CCRATIO.pdf")
p <- ggplot(results, aes(x=CC_RATIO, y=MEAN))+
  geom_point(size=3) +
  geom_smooth(method="lm", se=TRUE, fullrange=TRUE, level=0.95)
p <- p+geom_label(label=results$sites, nudge_y = 0.05, nudge_x=0.001)
p +theme_linedraw()
dev.off()