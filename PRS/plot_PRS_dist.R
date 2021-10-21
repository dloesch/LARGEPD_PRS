##read in phenotype file
pheno <- read.table("large_tb.merged.pheno.txt", header=TRUE)

#read in PRS file
PRS <- read.table("large_tb.merged.PRS.txt", header=TRUE)


#merge
data <- merge(pheno, PRS, by.x = "sample.id",by.y="ID", all.x=TRUE)

sites <- levels(data$SITE_STUDY)

#create groups for plotting
for(site in sites){
  
  n_cases <- nrow(data[data$PD_STATUS == 1 & data$SITE_STUDY == site,])
  n_controls <- nrow(data[data$PD_STATUS == 0 & data$SITE_STUDY == site,])
  
  if(n_cases== 0 | n_controls == 0){
    data$GROUP[data$SITE_STUDY == site] <- site
  }else{
    data$GROUP[data$PD_STATUS == 1 & data$SITE_STUDY == site] <- paste0(site, "_CASE")
    data$GROUP[data$PD_STATUS == 0 & data$SITE_STUDY == site] <- paste0(site, "_CONTROL")
  }
  
}

data$GROUP[data$PD_STATUS == 1 ] <- paste0(data$SITE_STUDY[data$PD_STATUS == 1], "_CASE")
data$GROUP[data$PD_STATUS == 0 ] <- paste0(data$SITE_STUDY[data$PD_STATUS == 0 ], "_CONTROL")
#create country label for fill

groups <- levels(data$GROUP)
means <- c()
for(g in groups){
  means <- c(means, mean(data$PRS_SCALED[data$GROUP == g]))
}

means <- as.data.frame(cbind(groups, means))
means$means <- as.numeric(as.character(means$means))
medians <- means[order(means$means),]

means <- means[order(means[,2]),]
data$GROUP <- factor(data$GROUP, means$groups)

for(i in 1:nrow(data)){
  data$COUNTRY[i] <- unlist(unname(strsplit(as.character(data$SITE_STUDY[i]), split="_")))[1]
}

#order by country
data <- data[order(data$COUNTRY),]
data$COUNTRY <- as.factor(data$COUNTRY)

#set up palette
colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


#plot 1
pdf("PRS_by_site.pdf")
p <- ggplot(data, aes(x=PRS_SCALED, y=GROUP, fill= COUNTRY)) + 
  geom_violin(trim=FALSE) + 
  stat_summary(fun.data=mean_sdl, mult=1, 
                 geom="pointrange", color="black") +
  theme_bw() +
  ggtitle("PRS: ordered by mean PRS")

p <- p + scale_fill_manual(values=colorBlindGrey8)
print(p)

#plot 2
data$GROUP <- factor(data$GROUP, levels = groups[order(groups)])
p <- ggplot(data, aes(x=PRS_SCALED, y=GROUP, fill= COUNTRY)) + 
  geom_violin(trim=FALSE) + 
  stat_summary(fun.data=mean_sdl, mult=1, 
               geom="pointrange", color="black") +
  theme_bw() +
  ggtitle("PRS: ordered by site")

p <- p + scale_fill_manual(values=colorBlindGrey8)
print(p)

#plot 3
p <- ggplot(data, aes(x=PRS_SCALED, y=COUNTRY, fill= COUNTRY)) + 
  geom_violin(trim=FALSE) + 
  stat_summary(fun.data=mean_sdl, mult=1, 
               geom="pointrange", color="black") +
  theme_bw() +
  ggtitle("PRS by country")+
  theme(legend.position = "none")   

p <- p + scale_fill_manual(values=colorBlindGrey8)


print(p)

dev.off()
