qsub -I -l ncpus=1,mem=5gb,cput=5:0:0 -l walltime=5:0:0 /bin/bash
module load intel intelmpi R
R


######################################################################################
#change input ONLY here
# ALL RUNS in CHERRY CREEK, FILES ARE IMPORTED AFTER PROCESSING

size <- matrix(,22,2)
size[,1]<- c(1:22)
size[,2]<- c(103,107,89,80,83,84,70,72,61,71,65,63,49,41,38,39,30,38,17,33,18,15)
size <- data.frame(size)
colnames(size)<-c("chr","chrlength")

for (k in 1:22){
  i=7
  chr=size$chr[k]
  length=size$chrlength[k]
  
  gs <- matrix(,length,5)
  
#########################################
#gold standard


for (i in 1:length){
  
  infile <- paste("/storage/nipm/kerimbae/pimass/output/revision/chr",chr,"/c",i,".mcmc.txt", sep = "")
  mcmc <- read.table(infile, head = T, stringsAsFactors = FALSE)
  
  gs[i,1] <- sum(mcmc$postc)
  gs[i,2] <- mean(mcmc$postc)
  gs[i,3] <- max(mcmc$postc)
  gs[i,4] <- min(mcmc$pos)
  gs[i,5] <- max(mcmc$pos)
}

colnames(gs) = c("sum","mean","max", "start","end")
write.table(gs, file =paste("/storage/nipm/kerimbae/pimass/output/routput/resultchr",chr,".txt",sep=""), row.names = F, col.names = T)

}


#locally 

chr = 2
mcmc <- read.table(paste0("C:/Users/Benazir/Desktop/revision/resultchr",chr,".txt"),header=T)

length=dim(mcmc)[1]
mcmc[1,]

mcmc[103,]


barplot(mcmc$sum, main=paste("Chromosome ", chr," Sum of Posterior Inclusion Probabilities"), xlab="Group number", names.arg= c(1:length))
abline(h = quantile(as.numeric(mcmc$sum), c(.9), na.rm=T), lty=2)
abline(h =quantile(as.numeric(mcmc$sum), c(.75), na.rm=T), col=c("red"))
legend(1, 1, c("90th percentile", "75th percentile"), col = c("black", "red"),lty = c(2, 1),lwd=c(1,1),bg="gray95")

#subset top 5% of the regions, then top 1% SNPs

# based on the true phenotype run routput/gschri.txt want to subset top 5% of the regions.
# we do it based on their rank and check if percentile is >=95
# percentiles for top based on 10K

### data about dataset cut by 100 SNPs each no overlap

size <- matrix(,22,2)
size[,1]<- c(1:22)
size[,2]<- c(103,107,89,80,83,84,70,72,61,71,65,63,49,41,38,39,30,38,17,33,18,15)
size <- data.frame(size)
colnames(size)<-c("chr","chrlength")

for (i in 1:22){
  
  chr=size$chr[i]
  chrlength=size$chrlength[i]
  results <- matrix(,chrlength,5)
  
  for (region in 1:chrlength){
    
    mcmc <- read.table(paste0("/storage/nipm/kerimbae/pimass/output/revision/chr",chr,"/c",region,".mcmc.txt"),header = T)
    summa = sum(mcmc$postc)
    start = mcmc$pos[1]
    end = mcmc$pos[dim(mcmc)[1]]
    
    results[region,] = rbind(chr, region,summa,start,end )
    
  }
  
  colnames(results) <-c("chr", "region","summa","start","end")
  
  results <- data.frame(results) 
  percentile <- ecdf(results$summa)(results$summa)
  ranks <- rank(-results$summa)
  grand <- cbind(results,percentile,ranks)
  
  write.table(grand, paste0("/storage/nipm/kerimbae/pimass/output/routput/summarychr",chr),col.names = T,row.names = T,sep = " ")
  
}


# combine 22 chrs in one data frame

grand <-  read.table(paste0("/storage/nipm/kerimbae/pimass/output/routput/summarychr1"),header=TRUE)

for (chr in 2:22){
  table <- read.table(paste0("/storage/nipm/kerimbae/pimass/output/routput/summarychr",chr),header=TRUE)
  grand <- rbind(grand,table)
  print(dim(grand))
  print(tail(grand))
  
}
write.table(grand,paste0("/storage/nipm/kerimbae/pimass/output/routput/grandsummary.csv"),sep = ",",row.names = TRUE,
            col.names = TRUE)

orderedGrand <- grand[order(-grand$summa),]
write.table(orderedGrand,paste0("/storage/nipm/kerimbae/pimass/output/routput/orderedGrandsummary.csv"),sep = ",",row.names = FALSE,
            col.names = TRUE)

### subsetting top 1% of SNPs 
grand <- subset(orderedGrand[1:64,])
length =dim(grand)[1]
results <- matrix(,length*10,7)
colnames(results) <-c("chr", "region","start","end","rsID","position","PIP")

for (i in 1:length) {
  
  chr = grand$chr[i]
  region = grand$region[i]
  
  mcmc <- read.table(paste0("/storage/nipm/kerimbae/pimass/output/revision/chr",chr,"/c",region,".mcmc.txt"),header = T)
  
  sorted <- mcmc[order(-mcmc$postc),]
  start= 10*(i-1)+1
  end=i*10
  results[start:end,1] <- chr
  results[start:end,2] <- region
  results[start:end,3] <- grand$start[i]
  results[start:end,4] <- grand$end[i]
  results[start:end,5] <- as.character(sorted$rs[1:10])
  results[start:end,6] <- sorted[1:10,3]
  results[start:end,7] <- sorted[1:10,4]
  
}

write.csv(results,"/storage/nipm/kerimbae/pimass/output/routput/table2overallwithPIP.csv")











