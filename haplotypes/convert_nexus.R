#script for converting haplotype results to nexus format

load("results/haplotypes/haplodata.RData")
popdata <- read.table("hap.populations.txt", header=FALSE)

#parameters
freq <- 0.005
prefix <- "rs356182"

data <- haplodata[[1]]

#convert data to counts
pops <- unique(popdata$V2)
for (p in pops){
    N <- nrow(popdata[popdata$V2 == p,])
    data[[p]] <- data[[p]]*N*2 
}

#fix Europe, one subject was dropped
N <- nrow(popdata[popdata$V2 == "EUR",])
data$EUR <- (data$EUR/(N*2))*(N-1)*2

#since this has a PD bias, will use filter based on 1KG
N <- sum(data$AFR)+sum(data$EUR)+sum(data$SAS) +sum(data$EAS) + sum(data$AMR)
data$KG <- (data$AFR + data$EUR + data$AMR + data$EAS + data$SAS)/N
data <- data[data$KG > freq,]

#data <- data[data$CASES > f | data$CONTROLS > f,] #for cases and controls

data$haplotypes <- as.character(data$haplotypes)
data$target <- ifelse(data$target == "REF", "G", "A")
data$hapID <- paste0(data$hapID, ":", data$target)

#file name
sink(paste0(prefix,".nex"))

#specify columns for annotations
superpops <- c("AFR",   "AMR", "EAS", "EUR", "SAS")
PD <- c("IPDGC_CASE",    "IPDGC_CONTROL", "LARGEPD_CASE", "LARGEPD_CONTROL")
PGP <- c("PGP_AMRIND", "PGP_MESTIZO")
#columns <-     c(superpops, PD, PGP)
columns <- c(superpops,PD)

#taxa block
cat("#NEXUS", "\n","\n" )
cat("BEGIN TAXA;","\n")
n <- nrow(data)
cat(paste0("DIMENSIONS NTAX=",n,";"),"\n","\n")
cat("TAXLABELS","\n")
for(i in 1:nrow(data)){
  cat(data$hapID[i],"\n")
}
cat(";","\n","\n")
cat("END;","\n","\n")

#Characters block
cat("BEGIN CHARACTERS;", "\n")
n <- nchar(data$haplotypes[1])
cat(paste0("DIMENSIONS NCHAR=",n,";"), "\n")
cat("FORMAT DATATYPE=DNA MISSING=? GAP=- ;", "\n")
cat("MATRIX", "\n" )
cat(write.table(data[c(1,2)], quote=FALSE, row.names=FALSE,
                col.names=FALSE, sep=" "))
cat(";", "\n" )
cat("END;", "\n", "\n")

#traits block, for adding colors
cat("BEGIN TRAITS;", "\n")
#specify which columns you want to include for traits
#n <- length(colnames(data)[columns])

n <- length(columns)
cat(paste0("  DIMENSIONS NTRAITS=",n,";"), "\n")
cat("  Format labels=yes missing=? separator=comma;", "\n")
#labels <- paste(colnames(data)[columns], collapse = " ")
labels <- paste(columns, collapse = " ")
cat(paste("  TraitLabels",labels),";", "\n")
cat("  Matrix", "\n")

#labels <- colnames(data)[columns]
labels <- columns

for(i in 1:nrow(data)){
  foo <- paste(data[labels][i,], collapse=",")
  cat(data$hapID[i],foo, "\n")
}
cat(";", "\n")
cat("END;", "\n")
sink()
