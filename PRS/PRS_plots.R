#!/usr/bin/Rscript

#script for calculating odds ratios by PRS, plotting odds ratios, and plotting AUC 

#packages
library(data.table)
library(tools)
library(boot)
library(cvTools)
library(pROC)
library(ggplot2)
library(caret)
library(dplyr)

#read in PRS file
prs.file <- "large_tb.merged.PRS.txt"
PRS <- read.table(prs.file, header=TRUE)


#read in pheno file
pheno.file <- "large_tb.merged.pheno.txt"


isCSV <- file_ext(pheno.file) == "csv"
if(isCSV){
  pheno <- read.csv(pheno.file, header=TRUE, comment.char = "#")
}else{
  pheno <- read.table(pheno.file, header=TRUE)
}


PRS_pheno <- merge(PRS, pheno, by.x="ID", by.y="sample.id", all.x = TRUE)


#prepare data
trait <- "PD_STATUS"


PRS_type <- "average"
prs <- ifelse(PRS_type == "raw", "PRS_RAW", ifelse(PRS_type == "average", "PRS_AVG", ifelse(PRS_type == "scaled", "PRS_SCALED", "PRS_AVG")))

#shuffle data
data <- PRS_pheno[sample(nrow(PRS_pheno)),]
data <- data[complete.cases(data),]

#formulas
f = as.formula(paste0(trait,"~",prs)) #PRS only

#Prediction via CV:
k <- 10
folds <- cvFolds(NROW(data), K=k) #original

data$PRS_only_pred <- rep(0,nrow(data))

for(i in 1:k){
  
  #original way of dividing data
  train <- data[folds$subsets[folds$which != i], ] 
  validation <- data[folds$subsets[folds$which == i], ] 
  
  new_PRS <- glm(f, data=train, family="binomial")
    
  #predictions
  PRS_only_pred <- predict.glm(new_PRS, newdata=validation, type="response")

  #add to data frame
  data[folds$subsets[folds$which == i], ]$PRS_only_pred <- PRS_only_pred
}

#generate confusion matrix
pred <- ifelse(data$PRS_only_pred > 0.5, 1, 0)


#AUC curves
prefix <- "large_tb.subset"
pdf(paste0(prefix,".ROC.pdf"))

#set up palettes
colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#roc of Peru
roc1 <- roc(response=data[[trait]][data$SITE == "Peru_Lima"], predictor = data$PRS_only_pred[data$SITE == "Peru_Lima"])
plot.roc(roc1, col=colorBlindGrey8[1], grid=TRUE)

#roc of large-PD
roc2 <- roc(response=data[[trait]][data$STUDY == "PD"], predictor = data$PRS_only_pred[data$STUDY == "PD"])
plot.roc(roc2,add=TRUE, col=colorBlindGrey8[7])

#roc of non-Peruvian large-PD
roc3 <- roc(response=data[[trait]][data$STUDY == "PD" & data$SITE != "Peru_Lima"], predictor = data$PRS_only_pred[data$STUDY == "PD" & data$SITE != "Peru_Lima"])
plot.roc(roc3,add=TRUE, col=colorBlindGrey8[3])

#roc of large-PD +TB
roc4 <- roc(response=data[[trait]], predictor = data$PRS_only_pred)
plot.roc(roc4,add=TRUE, col=colorBlindGrey8[6])


labels <- c(paste0("LARGE-PD: Peru-Lima (AUC:",signif(as.numeric(roc1$auc), digits = 3),")" ), 
            paste0("LARGE-PD (AUC:",signif(as.numeric(roc2$auc), digits = 3),")" ), 
            paste0("LARGE-PD: Non-Peruvian (AUC:",signif(as.numeric(roc3$auc), digits = 3),")" ), 
            paste0("LARGE-PD + Ext. Peruvian Controls (AUC:",signif(as.numeric(roc4$auc), digits = 3),")" ))
legend("bottomright", cex=0.8,lty=1, 
       col=c(colorBlindGrey8[1:3], colorBlindGrey8[6]), 
       legend=labels)

dev.off()
  
##odds ratios by PRS quintile

data$percentile <- as.factor(ntile(data$PRS_RAW, 100))
results <- as.data.frame(c("Q2", "Q3", "Q4", "Q5"))
data$quintile <- ifelse(data$percentile %in% 1:20, "Q1", 
                       ifelse(data$percentile %in% 21:40, "Q2", 
                              ifelse(data$percentile %in% 41:60, "Q3", 
                                     ifelse(data$percentile %in% 61:80, "Q4", "Q5"))))
tests <- c("Q2", "Q3", "Q4", "Q5")
ref <- "Q1"
bins <- "quintile"
colnames(results) <- bins
results$OR <- NA
results$lower <- NA
results$upper <- NA
results$ci <- NA

for(test in tests){
  Du <- nrow(data[data[[trait]] == 1 & data[[bins]] ==ref,])
  Hu <- nrow(data[data[[trait]] == 0 & data[[bins]] ==ref,])
  
  De <- nrow(data[data[[trait]] == 1 & data[[bins]] ==test,])
  He <- nrow(data[data[[trait]] == 0 & data[[bins]] ==test,])
  
  OR <- (De/He)/(Du/Hu)
  se <- 1.96*(sqrt(1/De + 1/He + 1/Du + 1/Hu))
  upper <- exp(log(OR) + se)
  lower <- exp(log(OR) - se)
  ci <- paste0(signif(lower, digits = 4),"-", signif(upper,digits=4))
  
  results$OR[results[[bins]] == test] <- OR
  results$lower[results[[bins]] == test] <- lower
  results$upper[results[[bins]] == test] <- upper
  results$ci[results[[bins]] == test] <- ci 
}

results$MODEL <- "LARGE-PD+EXT_CONTROLS"

temp <- results
data <- data[data$STUDY == "PD",]

data$percentile <- as.factor(ntile(data$PRS_RAW, 100))
results <- as.data.frame(c("Q2", "Q3", "Q4", "Q5"))
data$quintile <- ifelse(data$percentile %in% 1:20, "Q1", 
                        ifelse(data$percentile %in% 21:40, "Q2", 
                               ifelse(data$percentile %in% 41:60, "Q3", 
                                      ifelse(data$percentile %in% 61:80, "Q4", "Q5"))))
tests <- c("Q2", "Q3", "Q4", "Q5")
ref <- "Q1"
bins <- "quintile"
colnames(results) <- bins
results$OR <- NA
results$lower <- NA
results$upper <- NA
results$ci <- NA

for(test in tests){
  Du <- nrow(data[data[[trait]] == 1 & data[[bins]] ==ref,])
  Hu <- nrow(data[data[[trait]] == 0 & data[[bins]] ==ref,])
  
  De <- nrow(data[data[[trait]] == 1 & data[[bins]] ==test,])
  He <- nrow(data[data[[trait]] == 0 & data[[bins]] ==test,])
  
  OR <- (De/He)/(Du/Hu)
  se <- 1.96*(sqrt(1/De + 1/He + 1/Du + 1/Hu))
  upper <- exp(log(OR) + se)
  lower <- exp(log(OR) - se)
  ci <- paste0(signif(lower, digits = 4),"-", signif(upper,digits=4))
  
  results$OR[results[[bins]] == test] <- OR
  results$lower[results[[bins]] == test] <- lower
  results$upper[results[[bins]] == test] <- upper
  results$ci[results[[bins]] == test] <- ci 
}

results$MODEL <- "LARGE-PD"

results <- rbind(results, temp)

write.table(results, "results/PRS/PRS.OR.txt", sep='\t', quote=FALSE, row.names = FALSE, col.names = TRUE)

#prepare data
results$MODEL <- as.factor(results$MODEL)
results[[bins]] <- as.factor(results[[bins]])

g <- ggplot(results, aes(x =OR, y = results[[bins]], color=MODEL, fill=MODEL)) 

dodger = position_dodge(width = 0.3)
g <- g + geom_errorbarh(aes(xmax = upper, xmin = lower), size = .5, height = 
                          .2, color = "gray50", position = dodger) +
  geom_point(size = 3.5, position = dodger) +
  theme_bw() +
  coord_flip() +
  theme(panel.grid.minor = element_blank())  +
  ylab(toupper(bins)) +
  xlab("Odds Ratio (log10 scaled)") +
  ggtitle("PD Odds Ratios by PRS Quintile")

g <- g +scale_color_manual(values=c("orange", "dodgerblue"), name="Model", labels = c("LARGE-PD", "LARGE-PD + Ext. Controls"))

g <- g+guides(color=guide_legend("Model"), fill = FALSE)

g <- g+ scale_x_continuous(trans='log10')

pdf("large_tb.PRS_OR.pdf")
print(g)
dev.off()
