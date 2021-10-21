###Rscript for for comparing AUCs using Delong's method via pROC package
library(pROC)
library(cvTools)
#read in phenotypes
p1 <- data.table::fread("large.10PC.pheno.txt", data.table = FALSE)
p2 <- data.table::fread("large_tb.merged.10PC.pheno.txt", data.table = FALSE)
p3 <- data.table::fread("large.10PC.outliers.pheno.txt", data.table = FALSE)

#read in prs
prs1 <- data.table::fread("large.PRS.txt", data.table = FALSE)
prs2 <- data.table::fread("large.full_stats.PRS.txt", data.table=FALSE)
prs3 <- data.table::fread("large_tb.merged.PRS.txt", data.table=FALSE)

#merge
p1 <- merge(p1, prs1, by.x="sample.id", by.y="ID", all=FALSE, sort=FALSE)
colnames(prs2)[2] <- "PRS_FULL"
p1 <- merge(p1, prs2, by.x="sample.id", by.y="ID", all=FALSE, sort=FALSE)
p2 <- merge(p2, prs3, by.x="sample.id", by.y="ID", all=FALSE, sort=FALSE)
p3 <- merge(p3, prs1, by.x="sample.id", by.y="ID", all=FALSE, sort=FALSE)

#generate predictions
data <- p1[sample(nrow(p1)),]
data <- data[complete.cases(data),]

k <- 10
folds <- cvFolds(NROW(data), K=k) #original

data$pred <- rep(0,nrow(data))
trait <- "PD_STATUS"
prs <- "PRS_AVG"
f = as.formula(paste0(trait,"~",prs))

#GWAS significant, LARGE-PD only
for(i in 1:k){

  train <- data[folds$subsets[folds$which != i], ] 
  validation <- data[folds$subsets[folds$which == i], ] 
  
  fit <- glm(f, data=train, family="binomial")
    
  #predictions
  pred <- predict.glm(fit,newdata=validation, type="response")
    
  #add to data frame
  data[folds$subsets[folds$which == i], ]$pred <- pred
}


roc1 <- roc(response = data[[trait]], predictor = data$pred, ci=TRUE)

#Full sumstats, LARGE-PD only
prs <- "PRS_FULL"
f = as.formula(paste0(trait,"~",prs))
for(i in 1:k){
  
  train <- data[folds$subsets[folds$which != i], ] 
  validation <- data[folds$subsets[folds$which == i], ] 
  
  fit <- glm(f, data=train, family="binomial")
  
  #predictions
  pred <- predict.glm(fit,newdata=validation, type="response")
  
  #add to data frame
  data[folds$subsets[folds$which == i], ]$pred <- pred
}

roc2 <- roc(response = data[[trait]], predictor = data$pred, ci=TRUE)
test <- roc.test(roc1, roc2)

results <- data.frame(MODEL="FULL", PVAL=test$p.value)

#GWAS significant, LARGE-PD only, UNRELATED
unrel <- data.table::fread("large_tb.unrelated_toberemoved.txt", data.table=FALSE, header=FALSE)
data <- data[!data$sample.id %in% unrel$V2,]

folds <- cvFolds(NROW(data), K=k) 

data$pred <- rep(0,nrow(data))
trait <- "PD_STATUS"
prs <- "PRS_AVG"
f = as.formula(paste0(trait,"~",prs))


for(i in 1:k){
  
  train <- data[folds$subsets[folds$which != i], ] 
  validation <- data[folds$subsets[folds$which == i], ] 
  
  fit <- glm(f, data=train, family="binomial")
  
  #predictions
  pred <- predict.glm(fit,newdata=validation, type="response")
  
  #add to data frame
  data[folds$subsets[folds$which == i], ]$pred <- pred
}
roc3 <- roc(response = data[[trait]], predictor = data$pred, ci=TRUE)
test <- roc.test(roc1, roc3)

foo <- data.frame(MODEL="UNREL", PVAL=test$p.value)
results <- rbind(results, foo)

#GWAS significant, LARGE-PD only, downsampled
drop <- data.table::fread("drop.200_Peruvian_cases.txt", data.table = FALSE)
data <- p1[sample(nrow(p1)),]
data <- data[complete.cases(data),]
data <- data[!data$sample.id %in% drop$V1,]

folds <- cvFolds(NROW(data), K=k) 

data$pred <- rep(0,nrow(data))
trait <- "PD_STATUS"
prs <- "PRS_AVG"
f = as.formula(paste0(trait,"~",prs))


for(i in 1:k){
  
  train <- data[folds$subsets[folds$which != i], ] 
  validation <- data[folds$subsets[folds$which == i], ] 
  
  fit <- glm(f, data=train, family="binomial")
  
  #predictions
  pred <- predict.glm(fit,newdata=validation, type="response")
  
  #add to data frame
  data[folds$subsets[folds$which == i], ]$pred <- pred
}
roc4 <- roc(response = data[[trait]], predictor = data$pred, ci=TRUE)
test <- roc.test(roc1, roc4)

foo <- data.frame(MODEL="DOWN_SAMPLE", PVAL=test$p.value)
results <- rbind(results, foo)

##peru only
data <- p1[sample(nrow(p1)),]
data <- data[complete.cases(data),]
data <- data[data$SITE == "Peru",]

folds <- cvFolds(NROW(data), K=k) 

data$pred <- rep(0,nrow(data))
trait <- "PD_STATUS"
prs <- "PRS_AVG"
f = as.formula(paste0(trait,"~",prs))


for(i in 1:k){
  
  train <- data[folds$subsets[folds$which != i], ] 
  validation <- data[folds$subsets[folds$which == i], ] 
  
  fit <- glm(f, data=train, family="binomial")
  
  #predictions
  pred <- predict.glm(fit,newdata=validation, type="response")
  
  #add to data frame
  data[folds$subsets[folds$which == i], ]$pred <- pred
}
roc5 <- roc(response = data[[trait]], predictor = data$pred, ci=TRUE)
test <- roc.test(roc1, roc5)

foo <- data.frame(MODEL="PERU_ONLY", PVAL=test$p.value)
results <- rbind(results, foo)


#Peru excluded

data <- p1[sample(nrow(p1)),]
data <- data[complete.cases(data),]
data <- data[data$SITE != "Peru" & data$SITE != "Peru_Puno",]

folds <- cvFolds(NROW(data), K=k) 

data$pred <- rep(0,nrow(data))
trait <- "PD_STATUS"
prs <- "PRS_AVG"
f = as.formula(paste0(trait,"~",prs))


for(i in 1:k){
  
  train <- data[folds$subsets[folds$which != i], ] 
  validation <- data[folds$subsets[folds$which == i], ] 
  
  fit <- glm(f, data=train, family="binomial")
  
  #predictions
  pred <- predict.glm(fit,newdata=validation, type="response")
  
  #add to data frame
  data[folds$subsets[folds$which == i], ]$pred <- pred
}
roc6 <- roc(response = data[[trait]], predictor = data$pred, ci=TRUE)
test <- roc.test(roc1, roc6)

foo <- data.frame(MODEL="EX_PERU", PVAL=test$p.value)
results <- rbind(results, foo)


#large-PD + ext controls
data <- p2[sample(nrow(p2)),]
data <- data[complete.cases(data),]

folds <- cvFolds(NROW(data), K=k) #original

data$pred <- rep(0,nrow(data))

for(i in 1:k){
  
  train <- data[folds$subsets[folds$which != i], ] 
  validation <- data[folds$subsets[folds$which == i], ] 
  
  fit <- glm(f, data=train, family="binomial")
  
  #predictions
  pred <- predict.glm(fit,newdata=validation, type="response")
  
  #add to data frame
  data[folds$subsets[folds$which == i], ]$pred <- pred
}
roc7 <- roc(response = data[[trait]], predictor = data$pred, ci=TRUE)
test <- roc.test(roc1, roc7)
foo <- data.frame(MODEL="EXT_CONTROLS", PVAL=test$p.value)
results <- rbind(results, foo)

#large-PD + ext controls, unrelated
data <- data[!data$sample.id %in% unrel$V2,]
folds <- cvFolds(NROW(data), K=k) #original

data$pred <- rep(0,nrow(data))

for(i in 1:k){
  
  train <- data[folds$subsets[folds$which != i], ] 
  validation <- data[folds$subsets[folds$which == i], ] 
  
  fit <- glm(f, data=train, family="binomial")
  
  #predictions
  pred <- predict.glm(fit,newdata=validation, type="response")
  
  #add to data frame
  data[folds$subsets[folds$which == i], ]$pred <- pred
}
roc8 <- roc(response = data[[trait]], predictor = data$pred, ci=TRUE)
test <- roc.test(roc1, roc8)
foo <- data.frame(MODEL="EXT_CONTROLS_UNREL", PVAL=test$p.value)
results <- rbind(results, foo)

#large-PD + ext controls, peru only
data <- p2[sample(nrow(p2)),]
data <- data[complete.cases(data),]
data <- data[data$COUNTRY == "Peru",]

folds <- cvFolds(NROW(data), K=k) #original

data$pred <- rep(0,nrow(data))

for(i in 1:k){
  
  train <- data[folds$subsets[folds$which != i], ] 
  validation <- data[folds$subsets[folds$which == i], ] 
  
  fit <- glm(f, data=train, family="binomial")
  
  #predictions
  pred <- predict.glm(fit,newdata=validation, type="response")
  
  #add to data frame
  data[folds$subsets[folds$which == i], ]$pred <- pred
}
roc9 <- roc(response = data[[trait]], predictor = data$pred, ci=TRUE)
test <- roc.test(roc1, roc9)
foo <- data.frame(MODEL="EXT_CONTROLS_PERU", PVAL=test$p.value)
results <- rbind(results, foo)
write.table(results, "AUC_pval_delong.txt", sep='\t', quote=FALSE, row.names = FALSE, col.names = TRUE)


#add outliers model
#generate predictions
data <- p3[sample(nrow(p3)),]
data <- data[complete.cases(data),]

k <- 10
folds <- cvFolds(NROW(data), K=k) #original

data$pred <- rep(0,nrow(data))
trait <- "PD_STATUS"
prs <- "PRS_AVG"
f = as.formula(paste0(trait,"~",prs))


for(i in 1:k){
  
  train <- data[folds$subsets[folds$which != i], ] 
  validation <- data[folds$subsets[folds$which == i], ] 
  
  fit <- glm(f, data=train, family="binomial")
  
  #predictions
  pred <- predict.glm(fit,newdata=validation, type="response")
  
  #add to data frame
  data[folds$subsets[folds$which == i], ]$pred <- pred
}


roc10 <- roc(response = data[[trait]], predictor = data$pred, ci=TRUE)
test <- roc.test(roc1, roc10)
print(test)

#3rd degree relatives model
#generate predictions
data <- p1[sample(nrow(p1)),]
data <- data[complete.cases(data),]
data <- data[!data$sample.id %in% final,]
k <- 10
folds <- cvFolds(NROW(data), K=k) #original

data$pred <- rep(0,nrow(data))
trait <- "PD_STATUS"
prs <- "PRS_AVG"
f = as.formula(paste0(trait,"~",prs))


for(i in 1:k){
  
  train <- data[folds$subsets[folds$which != i], ] 
  validation <- data[folds$subsets[folds$which == i], ] 
  
  fit <- glm(f, data=train, family="binomial")
  
  #predictions
  pred <- predict.glm(fit,newdata=validation, type="response")
  
  #add to data frame
  data[folds$subsets[folds$which == i], ]$pred <- pred
}


roc11 <- roc(response = data[[trait]], predictor = data$pred, ci=TRUE)
test <- roc.test(roc1, roc11)
print(test)