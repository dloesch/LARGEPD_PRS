#this script is for resolving relative pairs using the output from king
#for input, this uses the king output that specifies the relative pairs, not the full kinship matrix
#for using the full kinship matrix, would need to process matrix first before using this script

kin.file <- "test.kin0"
kin <- data.table::fread(kin.file, data.table = FALSE)

#3rd degree
kin <- kin[kin$Kinship > 0.0442,]

#count instances of each subject
samples <- unique(c(kin$ID1, kin$ID2))
N <- c()
for(s in samples){
  N <- c(N, nrow(kin[kin$ID1 == s,])+nrow(kin[kin$ID2 == s,]))
}
samples <- data.frame(ID=samples, N=N, stringsAsFactors = FALSE)

temp_drop <- samples$ID[samples$N > 1]
drop <- c()

#filter out samples who occur more than once
for (i in sort(unique(samples$N), decreasing = TRUE)){
  drop <- c(drop,samples$ID[samples$N == i])
  kin <- kin[!kin$ID1 %in% drop,]
  kin <- kin[!kin$ID2 %in% drop,]
  
  if(length(temp_drop[temp_drop %in% unique(c(kin$ID1, kin$ID2))]) < 1){
   break
  }
}

#break ties usuing missingness from genotype file
temp.miss <- tempfile()
system2("plink", paste0("--bfile large.autosome.final --missing --out", temp.miss), stdout = FALSE)

#read in data, make two dataframes to merge join ID1 and ID2 fields
miss <- data.table::fread(paste0(temp.miss,".imiss"), data.table=FALSE)
id1 <- data.frame(ID=miss$IID, F1=miss$F_MISS)
id2 <- data.frame(ID=miss$IID, F2=miss$F_MISS)

#add missingness rate
kin <- merge(kin, id2, by.x="ID2", by.y="ID", all.x=TRUE, sort=FALSE)
kin <- merge(kin, id1, by.x="ID1", by.y="ID", all.x=TRUE, sort=FALSE)

#create second list with ties
drop2 <- ifelse(kin2$F1 > kin2$F2, kin2$ID1, kin2$ID2)

#create final list
final <- unique(c(drop, drop2))

#write out
write.table(final, "WORKSPACE/LARGE-PD/large.degree3.v1.txt", sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)