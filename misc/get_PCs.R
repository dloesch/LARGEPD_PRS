##script for running pcAiR and pc-Relate 

library(SeqArray)
library(SNPRelate)
library(gdsfmt)
library(tools)
library(GENESIS)
library(SeqVarTools)

gdsfile <- "./data/large_tb.pruned.pass2.gds"
prefix <- "large_tb"

#get ids for analysis
p <- read.table("large_tb.merged.pheno.txt", header=TRUE, stringsAsFactors=FALSE)
ids <- p$sample.id

#convert vcf to gds
fmt.import <- "GT"

#set reference genome
ref <- "hg38"


seqVCF2GDS(vcffile, gdsfile, fmt.import=fmt.import,reference=ref, 
           storage.option="LZMA_RA", parallel=TRUE, scenario = "general")

#convert gds to snp-gds
new_gds <- paste0(prefix,".pruned.snp.gds")
seqGDS2SNP(gdsfile,new_gds)

gdsfile <- new_gds
gds <- snpgdsOpen(gdsfile)

#run KING-robust algorithm#
gds <- snpgdsOpen(gdsfile)

ibd <- snpgdsIBDKING(gds,type="KING-robust", sample.id= ids,  autosome.only=FALSE, num.thread = countThreads())

outfile=paste(prefix,".ibd.RData",sep="")
save(ibd, file=outfile)

snpgdsClose(gds)

#run pcAIR
#parameters
sparse <- 0.01104854 #2^(-13/2), 5th degree
king <- kingToMatrix(ibd, thresh=sparse)
kinMat <- king
outfile <- paste(prefix,".pca.RData",sep="")
kin_thresh <- 0.04419417 # 2^(-9/2), 3rd degree
n_pcs <- 20
#min maf and missing rate
maf=0.01
miss=0.01

#gds file
gdsfile <- "./data/large_tb.pruned.pass2.gds"
gds <- seqOpen(gdsfile)

pca <- pcair(gds, kinobj=kinMat, divobj=king, num.cores=4,
             kin.thresh = kin_thresh, div.thresh= -kin_thresh, 
             maf=maf, missing.rate=miss, autosome.only=FALSE,
             sample.include=ids)

save(pca, file=outfile)

#pcrelate
seqData <- SeqVarData(gds)
var_block_size <- 1024
sample_block_size <- 10000
n_pcs <- 3

iterator <- SeqVarBlockIterator(seqData, variantBlock=sample_block_size)

pcrel <- pcrelate(iterator, sample.include=ids,
                  pcs=pca$vectors[,1:n_pcs],
                  training.set=pca$unrels,
                  maf.thresh=maf)

outfile <- paste(prefix,".pcrelate.RData", sep="")
save(pcrel, file=outfile)
seqClose(seqData)

##run 2nd iteration of pc-air
gdsfile <- "./data/large_tb.pruned.pass2.gds"
gds <- seqOpen(gdsfile)
kinMat <- pcrelateToMatrix(pcrel, scaleKin=2)
outfile <- paste(prefix,".pca2.RData",sep="")
pca <- pcair(gds, kinobj=kinMat, divobj=king, num.cores=4,
             kin.thresh = kin_thresh, div.thresh= -kin_thresh, 
             maf=maf, missing.rate=miss, autosome.only=FALSE,
             sample.include=ids)

save(pca, file=outfile)
seqClose(gds)
