#run test GWAS of large-PD + ext controls

library(Biobase)
library(GENESIS)
library(gdsfmt)
library(SeqVarTools)

#set prefix
prefix <- "large_tb"


#prepare phenotype file
pheno_file <- "large_tb.merged.10PC.pheno.txt"
pheno <- read.table(pheno_file, stringsAsFactors = FALSE, header=T, comment.char = "#")
outcome <- "PD_STATUS"
pheno$SEX <- as.factor(pheno$SEX)
covars <- c("AGE", "SEX", paste0("PC", 1:10))
sample.id <- pheno[,1]

#prepare annotated data frame
n <- ncol(pheno)
labels <- rep(NA,ncol(pheno))
metadata <- data.frame(labelDescriptions=labels)
annot <- AnnotatedDataFrame(pheno, metadata)

# kinship matrix or GRM
grm_file <- "results/PCA/large_tb.pcrelate.RData"
load(grm_file)
grm <- pcrelateToMatrix(pcrel, scaleKin=2)
grm_type <- "pcrelate"

#fit null model
family <- binomial
nullmod <- fitNullModel(annot, outcome=outcome, covars=covars,
                        cov.mat=grm, sample.id=sample.id,
                        family=family)
save(nullmod, file=paste(prefix,".null_model.RData", sep=""))

#test for association
#analysis characteristics
out_prefix = paste(prefix,".assoc_single", sep="")
test_type= "score"
variant_block_size= 1024
pass_only=FALSE #keep as false if using seqarray format
maf_threshold=0.01
mac_threshold=NA #MAC filtering not supported with snpgds format

#prepare genotype file
vcffile <- "data/large_tb.pruned.all.vcf.gz"
gdsfile <- "data/large_tb.pruned.all.gds"
seqVCF2GDS(vcffile, gdsfile, fmt.import="GT",reference="hg19", 
           storage.option="LZMA_RA", parallel=TRUE, scenario = "general")

#prepare seqdata file and iterator
gds <- seqOpen(gdsfile)
pData(annot)$sample.id <- as.character(pData(annot)$sample.id) #LRAGE-PD has numeric ids, change to character
#if genotype data has more subjects than phenotype data, use the next 3 lines
s <- as.data.frame(seqGetData(gds,"sample.id"),stringsAsFactors=FALSE)
colnames(s) <- "sample.id"
pData(annot) <- merge(s, pData(annot),by.x="sample.id",by.y="sample.id",all.x=TRUE,all.y=TRUE, sort=FALSE)
pData(annot) <- merge(s, pData(annot),by.x="sample.id",by.y="sample.id",all.x=TRUE,all.y=TRUE, sort=FALSE)
seqData <- SeqVarData(gds, sampleData=annot)


#filter
mac.min <- as.numeric(mac_threshold)
maf.min <- as.numeric(maf_threshold)

#the following code works for seqarray format gds files
if (!is.na(mac.min)) {
  filterByMAC(seqData, sample.id, mac.min=mac.min)
} 
#
if (!is.na(maf.min)) {
  filterByMAF(seqData, sample.id, maf.min=maf.min)
}

checkSelectedVariants(seqData)
block.size <- as.integer(variant_block_size)
iterator <- SeqVarBlockIterator(seqData, variantBlock=block.size)

#run association
test <- switch(tolower(test_type),
               score="Score",
               wald="Wald")
assoc <- assocTestSingle(iterator, nullmod, test=test, imputed=FALSE)

#save results
save(assoc, file=paste(prefix,".RData",sep=""))

#prepare qqplot
pvals <- assoc$Score.pval
#qqplot
obs <- sort(pvals)
obs <- obs[!is.na(obs)]
obs <- obs[is.finite(-log10(obs))]

exp <- c(1:length(obs))  / (length(1:length(obs))+1) 

x <- exp
y <- obs

thin <- FALSE
if (thin == TRUE){
  quant.subsample <- function(y, m=100, e=1) {
    # m: size of a systematic sample
    # e: number of extreme values at either end to use
    x <- sort(y)
    n <- length(x)
    quants <- (1 + sin(1:m / (m+1) * pi - pi/2))/2
    sort(c(x[1:e], quantile(x, probs=quants), x[(n+1-e):n]))
    # Returns m + 2*e sorted values from the EDF of y
  }
  n.x <- n.y <- length(obs)
  m <- .001 * max(n.x, n.y)
  e <- floor(0.0005 * max(n.x, n.y))
  x <- quant.subsample(exp, m, e)
  y <- quant.subsample(obs, m, e)
}


outfile <- paste0(prefix, ".qqplot.",trait,".pdf")
pdf(outfile)

plot(-log10(x), -log10(y),
     pch=".", cex=4,
     xlab="Expected p-values", ylab="Observed p-values", main="QQ Plot: GWAS")

abline(0,1,col="red",lwd=3, lty=3)

z = qnorm(pvals/2)
lambda = round(median(z^2) / 0.454, 3)
print(lambda)

legend('topleft', legend=paste0("GC Lambda: ",lambda))
dev.off()