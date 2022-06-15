##script for plotting PCA and PRS distribution
library(ggplot2)

#read in phenotype data
pheno <- read.table("large_tb.merged.pheno.txt",
                    header=TRUE, stringsAsFactors = FALSE)

#read in prs data
prs <- read.table("large_tb.merged.PRS.txt",
                  header=TRUE, stringsAsFactors = FALSE)

#merge
data <- merge(prs, pheno, by.x="ID", by.y="sample.id", all.x=FALSE, sort=FALSE)

#colorblind-friendly palette

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#additional colors
cbPalette2 <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499",
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
palette(cbPalette)

#PCA plot by quintile
data$GROUP <- as.factor(dplyr::ntile(data$PRS_AVG, 5))
pc <- ggplot(data[data$GROUP == 1 | data$GROUP == 5,], aes(x=PC1, y=PC2, colour=GROUP))+
  geom_point()
pc <- pc + scale_color_manual(values=cbPalette)
pc <- pc + labs(title="PCA by PRS Quntile")
#pc <- pc + theme(legend.position = "none")
print(pc)
pdf("./Desktop/PCA.dist_plots.PDF")

#PRS distrubtion by recruitment site
p1 <- ggplot(data, aes(PRS_SCALED, fill = COUNTRY)) +
  geom_density(alpha = 0.6) +
  xlim(min(data$PRS_SCALED)-1, max(data$PRS_SCALED)+1)+
  theme(axis.ticks.x=element_blank())+
  labs(fill="Country", x="Scaled Polygenic Risk Score", y=NULL, title="PRS distribution by country")
p1 <- p1+ scale_fill_manual(values=c("black", cbPalette2))
p1 <- p1+theme_bw()
print(p1)

#controls only
dat <- data[data$PD_STATUS == 0,]
p2 <- ggplot(dat, aes(PRS_SCALED, fill = COUNTRY)) +
  geom_density(alpha = 0.7) +
  xlim(min(dat$PRS_SCALED)-1, max(dat$PRS_SCALED)+1)+
  theme(axis.ticks.x=element_blank())+
  labs(fill="Country", x="Scaled Polygenic Risk Score", y=NULL, title="PRS distribution by Country: Controls")
p2 <- p2+ scale_fill_manual(values=cbPalette[c(1,3:8)])
p2 <- p2+theme_bw()
print(p2)


#PRS distribution in Peru
dat <- data[data$COUNTRY == "Peru",]
dat$GROUP <- ifelse(dat$PD_STATUS == 1, paste0(dat$SITE_STUDY, "_CASES"), dat$SITE_STUDY)
dat$GROUP <- ifelse(dat$GROUP == "Peru_PD", paste0(dat$GROUP, "_CONTROLS"), dat$GROUP)
dat$GROUP <- as.factor(dat$GROUP)

p3 <- ggplot(dat, aes(PRS_SCALED, fill = GROUP2)) +
  geom_density(alpha = 0.7) +
  xlim(min(dat$PRS_SCALED)-1, max(dat$PRS_SCALED)+1)+
  theme(axis.ticks.x=element_blank())+
  labs(fill="Peruvian Subset", x="Scaled Polygenic Risk Score", y=NULL, title="PRS distribution: Peruvians")
p3 <- p3+ scale_fill_manual(values=cbPalette2[9:12])
p3 <- p3+theme_bw()
print(p3)

#cluster by PC
dmat <- dist(scale(data[c(paste0("PC",1:2))]), method = "euclidean")
set.seed(240)

#hierarchical clustering, not used
#cl <- hclust(dmat, method="average")
#fit <- cutree(cl, k=6)
#data$CLUSTER <- as.factor(fit)

#kmeans clustering
k <- kmeans(dmat, iter.max = 100, nstart=10, centers=5)
data$CLUSTER <- as.factor(k$cluster)
plot(data$PC1, data$PC2, col=data$CLUSTER, pch=16)
legend("topleft", legend=paste0("Cluster",
                                levels(data$CLUSTER)), pch=16,col = levels(data$CLUSTER), cex=0.6)

#color PCA plot by cluser
pc1 <- ggplot(data, aes(x=PC1, y=PC2, colour=CLUSTER))+
  geom_point()
pc1 <- pc1 + scale_color_manual(values=cbPalette)
pc1 <- pc1 + labs(title="PCA: PC-derived clusters")
print(pc1)
pc1 <- pc1 + theme(legend.position = "none")


#PRS distrubtion by cluster
p4 <- ggplot(data, aes(PRS_SCALED, fill = CLUSTER)) +
  geom_density(alpha = 0.7) +
  xlim(min(data$PRS_SCALED)-1, max(data$PRS_SCALED)+1)+
  theme(axis.ticks.x=element_blank())+
  labs(fill="Cluster", x="Scaled Polygenic Risk Score", y=NULL, title="PRS distribution by Cluster")
p4 <- p4+ scale_fill_manual(values=cbPalette)
p4 <- p4+theme_bw()
print(p4)


#PRS distrubtion by cluster (controls only)
dat <- data[data$PD_STATUS == 0,]
p5 <- ggplot(dat, aes(PRS_SCALED, fill = CLUSTER)) +
  geom_density(alpha = 0.7) +
  xlim(min(dat$PRS_SCALED)-1, max(dat$PRS_SCALED)+1)+
  theme(axis.ticks.x=element_blank())+
  labs(fill="Cluster", x="Scaled Polygenic Risk Score", y=NULL, title="PRS distribution by cluster (clusters)")
p5 <- p5+ scale_fill_manual(values=cbPalette)
p5 <- p5+theme_bw()
print(p5)

#pca by cluster (controls only)
pca2 <- ggplot(dat, aes(x=PC1, y=PC2, colour=CLUSTER))+
  geom_point()
pca2 <- pca2 + scale_color_manual(values=cbPalette)
pca2 <- pca2 + labs(title="PCA: PC-derived clusters (controls)")
pca2 <- pca2 + theme_bw()
pca2 <- pca2 + theme(legend.position = "none")
print(pca2)

tiff(filename = "supplement.tiff", width=17, height=12, units="cm", res=300)
controls <- plot_grid(pca2, p5, nrow=1, labels = c("A.", "B."))
print(controls)
dev.off()

#PCA plot by Country
p <- ggplot(data, aes(x=PC1, y=PC2, colour=COUNTRY))+
  geom_point()
p <- p + scale_color_manual(values=c("black", cbPalette2))
p <- p + labs(title="PCA: Country")
print(p)
p <- p + theme(legend.position = "none")
dev.off()


#PCA plot: Peru
data$GROUP <- ifelse(data$COUNTRY == "Peru" & data$PD_STATUS == 1,
                     paste0(data$SITE_STUDY, "_CASE"),
                     ifelse(data$COUNTRY == "Peru" & data$PD_STATUS == 0,
                            paste0(data$SITE_STUDY, "_CONTROL"), "OTHER"))

data$GROUP <- as.factor(data$GROUP)

data$GROUP2 <- ifelse(data$GROUP == "Peru_PD_CASE", "CASES", 
                      ifelse(data$GROUP == "Peru_PD_CONTROL", "CONTROLS",
                             ifelse(data$GROUP == "Peru_Puno_PD_CONTROL", "PUNO",
                             ifelse(data$GROUP == "Peru_TB_CONTROL", "EXT.CONTROLS","OTHER"))))
data$GROUP2 <- as.factor(data$GROUP2)
pc2 <- ggplot(data[data$GROUP2 != "OTHER",], aes(x=PC1, y=PC2, colour=GROUP2))+
  geom_point()
#pc2 <- pc2 + scale_color_manual(values=c("white", cbPalette2[9:12]))
pc2 <- pc2 + scale_color_manual(values=cbPalette2[c(6,10:12)])
pc2 <- pc2 + labs(title="PCA: Peru")
print(pc2)
pc2 <- pc2 + theme(legend.position = "none")


#plot grid
plot_grid(p, p1,pc1,p4,
          labels = c("A.", "B.", "C.", "D."))

plot_grid(pc2, p3, labels=c("E.", "F."), nrow=2, ncol=2)

#CC ratio
results <- as.data.frame(unique(data$COUNTRY))
colnames(results) <-  "SITE"
for(s in results$SITE){
  results$N[results$SITE == s] <- nrow(data[data$COUNTRY == s,])
  results$CC[results$SITE == s] <- nrow(data[data$COUNTRY == s & data$PD_STATUS == 1,])/nrow(data[data$COUNTRY == s,])
  results$MEAN[results$SITE == s] <- mean(data$PRS_RAW[data$COUNTRY == s])
  results$SD[results$SITE == s] <- sd(data$PRS_RAW[data$COUNTRY == s])
}

print(results)
print(cor.test(results$CC[1:4]), results$MEAN[1:4])

#PCA plot by case-control status
data$PD_STATUS <- as.factor(data$PD_STATUS)
p <- ggplot(data, aes(x=PC1, y=PC2, colour=PD_STATUS))+
  geom_point()
p <- p + scale_color_manual(values=cbPalette)
p <- p + labs(title="PCA: Case-Control")
print(p)


plot_grid(p, p1,pc1,p4,
          labels = c("A.", "B.", "C.", "D."))

p3 <- ggplot(data[data$GROUP2 != "OTHER",], aes(PRS_SCALED, fill = GROUP2)) +
  geom_density(alpha = 0.7) +
  xlim(min(dat$PRS_SCALED)-1, max(dat$PRS_SCALED)+1)+
  theme(axis.ticks.x=element_blank())+
  labs(fill="Peruvian Subset", x="Scaled Polygenic Risk Score", y=NULL, title="PRS distribution: Peruvians")
p3 <- p3+ scale_fill_manual(values=cbPalette2[c(6,10:12)])
p3 <- p3+theme_bw()
print(p3)

p <- p + theme_classic() + theme(legend.position = "none")
pc1 <- pc1 + theme_classic() + theme(legend.position = "none")
pc2 <- pc2 + theme_classic()+ theme(legend.position = "none")

tiff(filename = "fig2.tiff", width=17, height=23.97, units="cm", res=300)
plot_grid(p, p1,pc1,p4, pc2, p3, labels=c("A.", "B.", "C.", "D.", "E.", "F."), nrow=3, ncol=2)
dev.off()
