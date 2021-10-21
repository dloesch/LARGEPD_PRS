#!/usr/bin/Rscript

library(argparser)
library(data.table)
library(tools)
library(boot)
library(cvTools)
library(DescTools)
library(pROC)
library(ggplot2)
library(caret)

#set argument parser
argp <- arg_parser("test PRS against phenotype data", hide.opts = TRUE)

#reguired arguments
argp <- add_argument(argp, "--PRS", help="PRS file", type = "character")
argp <- add_argument(argp, "--pheno", help="phenotype file", type = "character")
argp <- add_argument(argp, "--trait", help="phenotype of interest", type = "character")

#option arguments
argp <- add_argument(argp, "--prefix", help="output prefix", type = "character", default="output")

#optional file name columns
argp <- add_argument(argp, "--ID1", help="column name for PRS file IDs", type = "character", default = NA)
argp <- add_argument(argp, "--ID2", help="column name for phenotype file IDs",type = "character", default = NA)

#optional parameters
argp <- add_argument(argp, "--covars", help="list of covariates seperated by commas", type = "character", default = NA)
argp <- add_argument(argp, "--PRS_type", help="PRS type: raw, average, scaled", type = "character", default = "average")
argp <- add_argument(argp, "--prev", help="true prevalence of phenotype", default = NA)
argp <- add_argument(argp, "--emp", help="estimate empirical p-value", type= "logical", default = FALSE)
argp <- add_argument(argp, "--folds", help="number of folds for CV", type = "integer", default = 10)
argp <- add_argument(argp, "--trait_type", help="quantitative or binary", type="character", default="binary")
argp <- add_argument(argp, "--exclude", help="text file with subjct ids to exclude", type="character", default=NA)
argp <- add_argument(argp, "--include", help="text file with subjct ids to include", type="character", default=NA)

#read in phenotype file
#parse arguments
argsv <- parse_args(argp)


#read in PRS file
prs.file <- argsv$PRS
if(is.na(prs.file)){
  stop("error: please provide PRS file")
}else{
  PRS <- read.table(prs.file, header=TRUE)
}


id1 <- argsv$ID1
if(is.na(id1)){
  id1 <- colnames(PRS)[colnames(PRS) %in% c("ID", "id", "subjects", "SUBJECTS", "subject.id", "sample.id")]
  if(length(id1) != 1){
    stop("error: unable to infer PRS file ID column name, please provide")
  }
}


#read in pheno file
pheno.file <- argsv$pheno
if(is.na(pheno.file)){
  stop("error: please provide phenotype file")
}

isCSV <- file_ext(pheno.file) == "csv"
if(isCSV){
  pheno <- read.csv(pheno.file, header=TRUE, comment.char = "#")
}else{
  #pheno <- read.table(pheno.file, header=TRUE)
  pheno <- data.table::fread(pheno.file, data.table = FALSE)
}

id2 <- argsv$ID2
if(is.na(id2)){
  id2 <- colnames(pheno)[colnames(pheno) %in% c("ID", "id", "subjects", "SUBJECTS", "subject.id", "sample.id")]
  if(length(id2) != 1){
    stop("unable to infer phenotype file ID column name; please provide")
  }
}

PRS_pheno <- merge(PRS, pheno, by.x=id1, by.y=id2, all.x = TRUE)


#prepare data
trait <- argsv$trait
if(is.na(trait)){
  stop("error: please specify trait")
}

if(is.na(argsv$covars)){
  stop("error: please specify covariates to estimate variance")
}else{
  covariates <- unlist(unname(strsplit(argsv$covars, split=",")))
}

PRS_type <- argsv$PRS_type
prs <- ifelse(PRS_type == "raw", "PRS_RAW", ifelse(PRS_type == "average", "PRS_AVG", ifelse(PRS_type == "scaled", "PRS_SCALED", "PRS_AVG")))

#shuffle data
data <- PRS_pheno[sample(nrow(PRS_pheno)),]
data <- data[c(id1, trait, covariates, prs)]
data <- data[complete.cases(data),]

#drop subjects, if requested:
if(!is.na(argsv$exclude)){
  exclude.file <- argsv$exclude
  exclude <- read.table(exclude.file, header=FALSE, stringsAsFactors = FALSE)
  data <- data[!(data[[id1]] %in% exclude$V1),]
}

#include subjects, if requested:
if(!is.na(argsv$include)){
  include.file <- argsv$include
  include <- read.table(include.file, header=FALSE, stringsAsFactors = FALSE)
  data <- data[data[[id1]] %in% include$V1,]
}
#print number of subjects
print(paste("Testing", nrow(data), "subjects"))


#formulas
f = as.formula(paste0(trait,"~", paste(covariates, collapse = "+"),"+",prs)) #full
f2 = as.formula(paste0(trait,"~", paste(covariates, collapse = "+"))) #base
f3 = as.formula(paste0(trait,"~",prs)) #PRS only

trait_type <- argsv$trait_type

#fit 2 models
if(trait_type == "binary"){
  fit_full <-  glm(f, data = data, family = "binomial")
  fit_base <- glm(f2, data = data, family = "binomial")
  
}else{
  fit_full <-  lm(f, data = data)
  fit_base <- lm(f2, data = data)
  
}

#to get PRS p-value
print("estimating p-values and variance explained")
n <- nrow(summary(fit_full)$coefficients)
p_full <- summary(fit_full)$coefficients[n,4]


#variance explained on case control scale
#full model
if(trait_type == "binary"){
  var_full=as.numeric(PseudoR2(fit_full, which = "Nagelkerke") - PseudoR2(fit_base, which = "Nagelkerke"))
}else{
  var_full=as.numeric(summary(fit_full)$r.squared - as.numeric(fit_base)$r.squared)
}


#variance on liability scale
prev <- argsv$prev
if(!is.na(prev)){
  K <- as.numeric(prev) #prevalence 
  T=qnorm(1-K)
  z <- dnorm(T)
  P=nrow(data[data[[trait]] ==1,])/nrow(data)
  C=(K*(1-K)/z^2)*(K*(1-K))/(P*(1-P))
  full_l <- var_full*C
}else{
  full_l <- NA
}


#get empirical p-value
if(argsv$emp == TRUE){
  print("estimating empirical p-values")
  temp <- data
  
  indict <- 0
  
  for(i in 1:10000){
    temp[[trait]] <- sample(temp[[trait]], nrow(temp))
    
    if(trait_type == "binary"){
      fit1 <- glm(f, data = temp, family = "binomial")
    }else{
      fit1 <- lm(f, data = temp)
    }
    
    
    p1 <- summary(fit_full)$coefficients[n,4]
    
    if(p1 < p_full){
      indict <- indict +1
    }
    
  }
  
  p_full_emp <- (indict[1] +1)/(10000 +1) 
}else{
  p_full_emp <- NA
}


#CV error
print("estimating CV error")
if(trait_type == "binary"){
  cost <- function(r, pi = 0) mean(abs(r-pi) > 0.5)
  
  cv_error_10 <- cv.glm(glmfit=fit_full, cost=cost, data=data, K=10)$delta
}else{
 
  cv_error_10 <- cv.glm(glmfit=fit_full, data=data, K=10)$delta
}


#Prediction via CV:
print("predicting trait via 10-fold CV")

k <- argsv$folds
folds <- cvFolds(NROW(data), K=k) #original


data$PRS_pred <- rep(0,nrow(data))
data$base_pred <- rep(0,nrow(data))
data$PRS_only_pred <- rep(0,nrow(data))


for(i in 1:k){
  
  #original way of dividing data
  train <- data[folds$subsets[folds$which != i], ] 
  validation <- data[folds$subsets[folds$which == i], ] 
  
  
  if(trait_type == "binary"){
    new_fit <- glm(f,data=train, family="binomial") 
    new_base <- glm(f2, data=train, family="binomial")
    new_PRS <- glm(f3, data=train, family="binomial")
    
    #predictions
    PRS_pred <- predict.glm(new_fit,newdata=validation, type="response")
    base_pred <- predict.glm(new_base,newdata=validation, type="response")
    PRS_only_pred <- predict.glm(new_PRS, newdata=validation, type="response")
    
   
  }else{
    
    new_fit <- lm(f,data=train) 
    new_base <- lm(f2, data=train)
    new_PRS <- lm(f3, data=train)
    
    #predictions
    PRS_pred <- predict.lm(new_fit,newdata=validation)
    base_pred <- predict.lm(new_base,newdata=validation)
    PRS_only_pred <- predict.lm(new_PRS, newdata=validation)
    
   
  }

  
  #add to data frame
  data[folds$subsets[folds$which == i], ]$PRS_pred <- PRS_pred
  data[folds$subsets[folds$which == i], ]$base_pred <- base_pred
  data[folds$subsets[folds$which == i], ]$PRS_only_pred <- PRS_only_pred
  
 
}

#generate confusion matrix
pred <- ifelse(data$PRS_only_pred > 0.5, 1, 0)

cm <- confusionMatrix(table(pred,data$PD_STATUS))
acc <- unname(cm$overall[1])
bal_acc <- unname(cm$byClass[11])
sens <- unname(cm$byClass[1])
spec <- unname(cm$byClass[2])

#compile model stats
labels <- c("PVALUE", "EMP_P", "CV_ERROR", "R2", "R2_LIABILITY", "AUC", "AUC_CI", "ACC", "BAL_ACC", "SENS", "SPEC")
values <- c(p_full, p_full_emp, cv_error_10[2], var_full, full_l, NA, NA, acc, bal_acc, sens, spec)
model_stats <- cbind(labels, values)

prefix <- argsv$prefix
if(trait_type == "quantitative"){
  labels <- c(labels, "MSE", "RMSE")
  data$RES <- data$PRS_only_pred - data[[trait]]
  data$RES <- data$RES^2
  mse <- sum(data$RES)/nrow(data)
  rmse <- sqrt(mse)
  values <- c(values, mse, rmse)
  model_stats <- cbind(labels, values)
  write.table(model_stats, paste0(prefix, ".model_stats.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}else{
  print("generating AUC and ROC plots")
  
  #AUC curves
  pdf(paste0(prefix,".ROC.pdf"))
  roc1 <- roc(response = data[[trait]], predictor = data$PRS_pred, ci=TRUE)
  plot.roc(roc1, print.auc.x = 1, print.auc.y = 1,print.auc.col = "dodgerblue", col="dodgerblue")
  roc2 <- roc(response=data[[trait]], predictor = data$base_pred, ci=TRUE)
  plot.roc(roc2, add = TRUE,print.auc.x=1,print.auc.y = 0.95)
  roc3 <- roc(response=data[[trait]], predictor = data$PRS_only_pred, ci=TRUE)
  plot.roc(roc3, add = TRUE, print.auc.x = 1, print.auc.y = 0.9, col= "tomato3", print.auc.col = "tomato3")
  legend("bottomright", lty=1, col=c("dodgerblue", "black", "tomato3"), legend=c("full model", "base model", "PRS only"))
  dev.off()
  
  #AUC values
  ci1 <- paste0(signif(roc3$ci[1], digits = 3),"-",signif(roc3$ci[3], digits=3))
  model_stats[6:7,2] <- c(roc3$auc, ci1)
  write.table(model_stats, paste0(prefix, ".model_stats.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
              
  #plot case-control separation:
  pdf(paste0(prefix,".distribution.pdf"))
  data[[trait]] <- ifelse(data[[trait]] == 1, "CASE", "CONTROL")
  data[[trait]] <- as.factor(data[[trait]])
  
  if(PRS_type != "scaled"){
    data$PRS <- scale(data[[prs]])
  }else{
    data$PRS <- data[[prs]]
  }
  print("plotting PRS distribution")
  p1 <- ggplot(data, aes(PRS, fill = data[[trait]])) +
    geom_density(alpha = 0.8) +
    xlim(min(data$PRS)-1, max(data$PRS)+1)+
    theme(axis.ticks.x=element_blank())+
    labs(fill=paste(trait), x="Scaled Polygenic Risk Score", y=NULL, title="PRS distribution: cases versus controls")
  p1 <- p1+ scale_fill_manual(values=c("dodgerblue", "orange", "white"))
  p1 <- p1+theme_bw()
  print(p1)
  dev.off()
}





