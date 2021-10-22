##script for survival analysis of PD AAO using survival package in R
##Kaplan-meier curves and cox regression stratified by PRS quintile
#note: script was modified to remove external controls. Results remain the same but makes methods more straightforward
#load packages
library(survival)
library(dplyr)
library(survminer)

#read in PRS data
#prs <- read.table("large_tb.merged.PRS.txt", header=TRUE, stringsAsFactors = FALSE) #external controls
prs <- read.table("large.PRS.txt", header=TRUE, stringsAsFactors = FALSE)

# read in phenotpe data
#pheno <- read.table("large_tb.merged.10PC.pheno.txt", header=TRUE, stringsAsFactors = FALSE) #external controls
pheno <- read.table("large.10PC.pheno.txt", header=TRUE, stringsAsFactors = FALSE)

#read in AAO data
aao <- read.csv("large-PD.AGE_ONSET.csv", header=TRUE)


#create file merging pheno, aao, and PRS data
p <- merge(prs, pheno, by.x="ID", by.y="sample.id", all.x = TRUE)
p <- merge(p, aao, by.x="ID", by.y="ID", all.x=TRUE)
p$AGE_ONSET <- ifelse(is.na(p$AGE_ONSET), p$AGE, p$AGE_ONSET)

#filter data for unrelated individuals
rels <- read.table("large_tb.unrelated_toberemoved.txt", header=FALSE)
p <- p[!p$ID %in% rels$V1,]

#remove outliers by AAO (drop if AA0 is < 18)
p <- p[p$AGE_ONSET > 18,]
p <- p[complete.cases(p),]

#stratify PRS by bin (Quintile, Quartile, Decile)
#will use quintile for manuscript as that matches papers I have recently read
p$percentile <- ntile(p$PRS_RAW, 100)
p$QUARTILE <- ifelse(p$percentile > 75, "Q4", 
                     ifelse( p$percentile > 50 & p$percentile <= 75, "Q3", 
                             ifelse( p$percentile > 25 & p$percentile <=50, "Q2", "Q1")))


p$QUINTILE <- ifelse(p$percentile %in% 1:20, "Q1", 
                     ifelse(p$percentile %in% 21:40, "Q2", 
                            ifelse(p$percentile %in% 41:60, "Q3", 
                                   ifelse(p$percentile %in% 61:80, "Q4", "Q5"))))

p$DECILE <- ifelse(p$percentile %in% 1:10, "P1", 
                   ifelse(p$percentile %in% 11:20, "P2", 
                          ifelse(p$percentile %in% 21:30, "P3", 
                                 ifelse(p$percentile %in% 31:40, "P4",
                                        ifelse(p$percentile %in% 41:50, "P5",
                                               ifelse(p$percentile %in% 51:60, "P6",
                                                      ifelse(p$percentile %in% 61:70, "P7", 
                                                             ifelse(p$percentile %in% 71:80, "P8",
                                                                    ifelse(p$percentile %in% 81:90, "P9","P10")))))))))

#set pallette for plots
colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#set group
group <- "QUINTILE"
p$GROUP <- as.factor(p[[group]])

#cases only
p1 <- p[p$PD_STATUS ==1,]
p1$surv <- Surv(time=p1$AGE_ONSET, p1$PD_STATUS)

#just plot top and bottom bins  
if(group == "QUARTILE"){
  p1$GROUP <- as.factor(ifelse(p1$GROUP == "Q1", "BOTTOM", ifelse(p1$GROUP == "Q4", "TOP", NA)))
}

if(group == "QUINTILE"){
  p1$GROUP <- as.factor(ifelse(p1$GROUP == "Q1", "BOTTOM", ifelse(p1$GROUP == "Q5", "TOP", NA)))
}

if(group == "DECILE"){
  p1$GROUP <- as.factor(ifelse(p1$GROUP == "P1", "BOTTOM", ifelse(p1$GROUP == "P10", "TOP", NA)))
}

#Kaplan-meier curver
kms <- survfit(surv ~ GROUP, data=p1[!is.na(p1$GROUP),])
g <- ggsurvplot(kms, conf.int = TRUE, surv.median.line = c("hv"), data=p1[!is.na(p1$GROUP),], pval = TRUE, risk.table = FALSE, 
           title=paste0("AAO by PRS ", group, ": PD Cases Only"),
           xlab = "Age at PD Diagnosis", ylab="Probability", palette=colorBlindGrey8[2:8], pval.method = TRUE, break.x.by=5, xlim=c(0,85))

pdf("large-PD.KM.quintile.pdf")
print(g)
dev.off()


#now all data (cases + controls)
p$GROUP <- as.factor(p[[group]])
p1 <- p
p1$surv <- Surv(time=p1$AGE_ONSET, p1$PD_STATUS)
#just plot top and bottom
if(group == "QUARTILE"){
  p1$GROUP <- as.factor(ifelse(p1$GROUP == "Q1", "BOTTOM", ifelse(p1$GROUP == "Q4", "TOP", NA)))
}

if(group == "QUINTILE"){
  p1$GROUP <- as.factor(ifelse(p1$GROUP == "Q1", "BOTTOM", ifelse(p1$GROUP == "Q5", "TOP", NA)))
}

if(group == "DECILE"){
  p1$GROUP <- as.factor(ifelse(p1$GROUP == "P1", "BOTTOM", ifelse(p1$GROUP == "P10", "TOP", NA)))
}
kms <- survfit(surv ~ GROUP, data=p1[!is.na(p1$GROUP),])
g <- ggsurvplot(kms, conf.int = TRUE, surv.median.line = c("hv"), data=p1[!is.na(p1$GROUP),], pval = TRUE, risk.table = FALSE, 
                title=paste0("AAO by PRS ", group, ": All Subjects"),
                xlab = "Age at PD Diagnosis", ylab="Probability", palette=colorBlindGrey8[2:8], pval.method = TRUE, break.x.by=5, xlim=c(0,85))

#pdf("large_tb.KM.quintile.pdf")
pdf("large.ALL.KM.quintile.pdf")
print(g)
dev.off()

##cox regression
#PD-cases only
group <- "QUINTILE"
p$GROUP <- as.factor(p[[group]])
p1 <- p[p$PD_STATUS ==1,]
p1$surv <- Surv(time=p1$AGE_ONSET, p1$PD_STATUS)
#coxreg1 <- coxph(surv~GROUP+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+SEX+SITE_STUDY, data=p1)
coxreg1 <- coxph(surv~GROUP+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+SEX+SITE, data=p1)
results <- as.data.frame(summary(coxreg1)$coefficients[1:4,])

#add confidence intervals
ci <- summary(coxreg1)$conf.int[1:4,]
results$lower <- unlist(unname(ci[,3]))
results$upper <- unlist(unname(ci[,4]))
results$ci <- paste(signif(results$lower,3), "-", signif(results$upper,3))

#non-peru sites
#coxreg1 <- coxph(surv~GROUP+PC1+PC2+PC3+PC4+PC5+SEX+SITE_STUDY, data=p1[p1$SITE != "Peru_Lima",])
coxreg1 <- coxph(surv~GROUP+PC1+PC2+PC3+PC4+PC5+SEX+SITE, data=p1[p1$SITE != "Peru_Lima",])
results_np <- as.data.frame(summary(coxreg1)$coefficients[1:4,])

ci <- summary(coxreg1)$conf.int[1:4,]
results_np$lower <- unlist(unname(ci[,3]))
results_np$upper <- unlist(unname(ci[,4]))
results_np$ci <- paste(signif(results$lower,3), "-", signif(results$upper,3))

#peru sites
#coxreg1 <- coxph(surv~GROUP+PC1+PC2+PC3+PC4+PC5+SEX, data=p1[p1$SITE == "Peru_Lima",])
coxreg1 <- coxph(surv~GROUP+PC1+PC2+PC3+PC4+PC5+SEX, data=p1[p1$SITE == "Peru",])
results_p <- as.data.frame(summary(coxreg1)$coefficients[1:4,])

ci <- summary(coxreg1)$conf.int[1:4,]
results_p$lower <- unlist(unname(ci[,3]))
results_p$upper <- unlist(unname(ci[,4]))
results_p$ci <- paste(signif(results$lower,3), "-", signif(results$upper,3))

###all subjects
p1 <- p
p1$surv <- Surv(time=p1$AGE_ONSET, p1$PD_STATUS)
#coxreg1 <- coxph(surv~GROUP+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+SEX+SITE_STUDY, data=p1)
coxreg1 <- coxph(surv~GROUP+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+SEX+SITE, data=p1)
results2 <- as.data.frame(summary(coxreg1)$coefficients[1:4,])

ci <- summary(coxreg1)$conf.int[1:4,]
results2$lower <- unlist(unname(ci[,3]))
results2$upper <- unlist(unname(ci[,4]))
results2$ci <- paste(signif(results2$lower,3), "-", signif(results2$upper,3))

#save
write.table(results, "large-PD.coxreg_results.cases.10PC.txt", sep='\t', quote=FALSE, row.names = TRUE, col.names = TRUE)
#write.table(results2, "large-PD_tb.coxreg_results.all_subjects.10PC.txt", sep='\t', quote=FALSE, row.names = TRUE, col.names = TRUE)
write.table(results2, "large-PD.coxreg_results.all_subjects.10PC.txt", sep='\t', quote=FALSE, row.names = TRUE, col.names = TRUE)
write.table(results_p, "large-PD.coxreg_results.peru.txt", sep='\t', quote=FALSE, row.names = TRUE, col.names = TRUE)
write.table(results_np, "large-PD.coxreg_results.other.txt", sep='\t', quote=FALSE, row.names = TRUE, col.names = TRUE)
