###script for plotting haplotypes###
library(ggplot2)
library(reshape2)
source("scripts/hap_functions.R")
#load data
load("results/haplotypes/haplodata.RData")
load("results/haplotypes/haplotypes.RData")

#plot
d <- haplodata[[1]]
data <- haplodata[[2]]
labels <- haplodata[[3]]


#reformat data by haplotype
ref <- unique(labels$topREF)
alt <- unique(labels$topALT)
pops <- as.character(data$pops)

for(h in c(ref,alt)){
  for(p in pops){
    data[[h]][data$pops == p] <- d[[p]][d$hapID ==h]
  }
}

#re-calculate other ref and alt columns:
data$otherREF <- (data$topREF_freq + data$otherREF) - (data$hap1 + data$hap9)
data$otherALT <- (data$topALT_freq + data$otherALT) - (data$hap2 + data$hap52 + data$hap11)
data$topALT_freq <- NULL
data$topREF_freq <- NULL


#format for plotting
tbl <- t(as.matrix(data[2:ncol(data)]))
tbl <- as.data.frame(tbl)
colnames(tbl) <- data$pops


#tbl$row <- as.factor(rownames(tbl))
tbl$row <- rownames(tbl)
dat2 <- melt(tbl, id.vars = "row")
dat2$row <- factor(dat2$row, levels=c(ref, "otherREF", alt, "otherALT"))


for(i in 1:nrow(dat2)){
  dat2$label[i] <- ifelse(dat2$row[i] == "topREF_freq", labels$topREF[labels$pops == dat2$variable[i]], ifelse(dat2$row[i] == "topALT_freq",
                                                                                                               labels$topALT[labels$pops == dat2$variable[i]], as.character(dat2$row[i])))
}

#dat2$label <- factor(dat2$label, levels=(c(unique(labels$topREF), "otherREF", unique(labels$topALT), "otherALT")))
dat2$label <- dat2$row

key <- c(paste0("REF haplotype: ", ref, "\n"), "other REF haplotypes\n", paste0("ALT haplotype: ", alt, "\n"), "other ALT haplotypes \n")

plot_colors <- c("firebrick4", "firebrick2", "slategray4","steelblue1", "steelblue3", "steelblue4","slategray3")

snp <- "rs356182"
key <- gsub("REF", "G allele \n", key)
key <- gsub("ALT", "A allele \n", key)

pdf("rs356182_haplotypes.PDF")

dat2$KEY <- dat2$label
ggplot(dat2, aes(x=variable, y=value, fill=KEY)) +
  geom_bar(stat="identity") +
  xlab("Population") +
  ylab(paste0(snp, " haplotypes")) +
  theme_classic() +
  scale_fill_manual(values = plot_colors, labels=key) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle(paste("Global", snp, "haplotypes by population"))



#matrix of shared alleles
haplomat <- get_haplomat(haplodata)

melted_cormat <- prep_haplomat(haplomat)

ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "white", high = "firebrick2", mid = "lightblue", 
                       midpoint = min(melted_cormat$value), limit = c(0,1), space = "Lab", 
                       name="Similarity") +
  theme_minimal()+ # minimal theme
  geom_text(aes(Var2, Var1, label = signif(value, digits = 2)), color = "black", size = 4) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 9, hjust = 1),
        axis.title.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.justification = c(1, 0),
        legend.position = c(0.4, 0.6),
        legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))+
  coord_fixed()+ 
  ggtitle(paste("Proportion of shared variant alleles between global", snp, "haplotypes"))

dev.off()

#generate summary table
haplotypes <- unique(c(haplodata[[3]]$topREF, haplodata[[3]]$topALT))
data <- haplodata[[1]][haplodata[[1]]$hapID %in% haplotypes,]
data$haplotypes <- NULL

write.table(data, "rs356182_haplotypes.summary.txt", sep='\t', quote=FALSE, row.names = FALSE, col.names = TRUE)
