####PRS distribution and Nalls et al. allele frequencies in 1000 Genomes####
library(ggplot2)

#pop labels
kg <- read.table("integrated_call_samples_v3.20130502.ALL.panel", header=TRUE, stringsAsFactors=FALSE)

#allele frequencies
freq <- read.table("results/1KG.freq_PD_snps.txt", header=TRUE, stringsAsFactors =FALSE)

#nalls sumstats for direction of effect
nalls <- read.table("nalls.PD.sumstats.txt", header=TRUE, stringsAsFactors=FALSE)

#PRS
prs <- read.table("1KG.PRS.txt", header=TRUE, stringsAsFactors=FALSE)

#merge 1KG pop file and PRS file
prs <- merge(prs, kg, by.x="ID", by.y="sample", all.x=FALSE, all.y=TRUE)

#merge nalls sumstats and superpop frequencies
freq <- freq[c("snps", "EUR", "EAS", "AMR", "SAS", "AFR")]
freq <- merge(freq, nalls, by.x="snps", by.y="SNP", all.x=TRUE)

#wilcoxon for PRS by superpop
#set EUR as reference

pops <- unique(prs$super_pop)
data <- as.data.frame(pops)
colnames(data) <- "POP"

for(pop in pops){
  print(pop)
  data$MEAN[data$POP == pop] <- mean(prs$PRS_RAW[prs$super_pop == pop])
  data$SD[data$POP == pop] <- sd(prs$PRS_RAW[prs$super_pop == pop])
  test <- wilcox.test(prs$PRS_RAW[prs$super_pop == "EUR"], prs$PRS_RAW[prs$super_pop == pop], alternative = "two.sided")
  data$PVAL[data$POP == pop] <- test$p.value
}

write.table(data,"PRS.wilcoxon.txt", sep='\t', quote = FALSE, row.names = FALSE, col.names = TRUE)


#plot PRS by superpop
prs$super_pop <- as.factor(prs$super_pop)

colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

pdf("1KG.PD_PRS.pdf")
p <- ggplot(prs, aes(x=PRS_SCALED, y=super_pop, fill= super_pop)) + 
  geom_violin(trim=FALSE) + 
  stat_summary(fun.data=mean_sdl, mult=1, 
               geom="pointrange", color="black") +
  theme_bw() +
  ggtitle("PRS by 1000 Genomes Super-Population") +
  theme(legend.position = "none") 

p <- p + scale_fill_manual(values=colorBlindGrey8[2:length(colorBlindGrey8)])
print(p)
dev.off()

#plot allele frequency by super_pop
freq$EFFECT <- ifelse(freq$BETA > 0, "RISK", "PROTECTIVE")
freq$EFFECT <- factor(freq$EFFECT, levels=c("RISK", "PROTECTIVE"))
pdf("PD_AF.plots.pdf")
f <- ggplot(freq, aes(x=AFR, y=EUR, color=EFFECT))+
  geom_point(size=3)+
  geom_abline(intercept=0, slope=1)+
  theme_classic()
f <- f+scale_color_manual(values=colorBlindGrey8[c(2,6)])+
  ggtitle("PD Risk allele frequencies: AFR vs EUR")
print(f)

f <- ggplot(freq, aes(x=EAS, y=EUR, color=EFFECT))+
  geom_point(size=3)+
  geom_abline(intercept=0, slope=1)+
  theme_classic()
f <- f+scale_color_manual(values=colorBlindGrey8[c(2,6)])+
  ggtitle("PD Risk allele frequencies: EAS vs EUR ")
print(f)

f <- ggplot(freq, aes(x=AMR, y=EUR, color=EFFECT))+
  geom_point(size=3)+
  geom_abline(intercept=0, slope=1)+
  theme_classic()
f <- f+scale_color_manual(values=colorBlindGrey8[c(2,6)])+
  ggtitle("PD Risk allele frequencies: AMR vs EUR ")
print(f)

f <- ggplot(freq, aes(x=SAS, y=EUR, color=EFFECT))+
  geom_point(size=3)+
  geom_abline(intercept=0, slope=1)+
  theme_classic()
f <- f+scale_color_manual(values=colorBlindGrey8[c(2,6)])+
  ggtitle("PD Risk allele frequencies: SAS vs EUR ")
print(f)
dev.off()


#construct contingency tables
sink("PD_AF.chisqr.txt")
#col: higher AF EUR vs AFR, TRUE FALSE
#row: Risk: TRUE FALSE
afr <- matrix(nrow = 2, ncol = 2)
afr[1,1] <- nrow(freq[freq$EFFECT == "RISK" & freq$EUR > freq$AFR,])
afr[1,2] <- nrow(freq[freq$EFFECT == "RISK" & freq$EUR < freq$AFR,])
afr[2,1] <- nrow(freq[freq$EFFECT != "RISK" & freq$EUR > freq$AFR,])
afr[2,2] <- nrow(freq[freq$EFFECT != "RISK" & freq$EUR < freq$AFR,])
colnames(afr) <- c("EUR>AFR", "EUR<AFR")
rownames(afr) <- c("RISK", "PROTECTIVE")
print(afr)
chisq.test(afr)

#EAS
#col: higher AF EUR vs EAS, TRUE FALSE
#row: Risk: TRUE FALSE
eas <- matrix(nrow = 2, ncol = 2)
eas[1,1] <- nrow(freq[freq$EFFECT == "RISK" & freq$EUR > freq$EAS,])
eas[1,2] <- nrow(freq[freq$EFFECT == "RISK" & freq$EUR < freq$EAS,])
eas[2,1] <- nrow(freq[freq$EFFECT != "RISK" & freq$EUR > freq$EAS,])
eas[2,2] <- nrow(freq[freq$EFFECT != "RISK" & freq$EUR < freq$EAS,])
colnames(eas) <- c("EUR>EAS", "EUR<EAS")
rownames(eas) <- c("RISK", "PROTECTIVE")
print(eas)
chisq.test(eas)

#AMR
amr <- matrix(nrow = 2, ncol = 2)
amr[1,1] <- nrow(freq[freq$EFFECT == "RISK" & freq$EUR > freq$AMR,])
amr[1,2] <- nrow(freq[freq$EFFECT == "RISK" & freq$EUR < freq$AMR,])
amr[2,1] <- nrow(freq[freq$EFFECT != "RISK" & freq$EUR > freq$AMR,])
amr[2,2] <- nrow(freq[freq$EFFECT != "RISK" & freq$EUR < freq$AMR,])
colnames(amr) <- c("EUR>AMR", "EUR<AMR")
rownames(amr) <- c("RISK", "PROTECTIVE")
print(amr)
chisq.test(amr)
#SAS
sas <- matrix(nrow = 2, ncol = 2)
sas[1,1] <- nrow(freq[freq$EFFECT == "RISK" & freq$EUR > freq$SAS,])
sas[1,2] <- nrow(freq[freq$EFFECT == "RISK" & freq$EUR < freq$SAS,])
sas[2,1] <- nrow(freq[freq$EFFECT != "RISK" & freq$EUR > freq$SAS,])
sas[2,2] <- nrow(freq[freq$EFFECT != "RISK" & freq$EUR < freq$SAS,])
colnames(sas) <- c("EUR>SAS", "EUR<SAS")
rownames(sas) <- c("RISK", "PROTECTIVE")
print(sas)
chisq.test(sas)

sink()
