# Table 2 new version top 5% overall regions 
qsub -I -l ncpus=1,mem=10gb,cput=5:0:0 -l walltime=5:0:0 /bin/bash
module load intel intelmpi R
R

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
    
    mcmc <- read.table(paste0("/storage/nipm/kerimbae/sliding_window/chr",chr,"/c",region,".mcmc.txt"),header = T)
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
  
  write.table(grand, paste0("/storage/nipm/kerimbae/sliding_window/summarychr",chr),col.names = T,row.names = T,sep = " ")
  
}


# combine 22 chrs in one data frame

grand <-  read.table(paste0("/storage/nipm/kerimbae/sliding_window/summarychr1"),header=TRUE)

for (chr in 2:22){
  table <- read.table(paste0("/storage/nipm/kerimbae/sliding_window/summarychr",chr),header=TRUE)
  grand <- rbind(grand,table)
  print(dim(grand))
  print(tail(grand))
  
}
write.table(grand,paste0("/storage/nipm/kerimbae/sliding_window/grandsummary.csv"),sep = ",",row.names = TRUE,
            col.names = TRUE)

# sort table in excel don't skip this step


### subsetting top 1% of SNPs 
grand <- read.csv("/storage/nipm/kerimbae/sliding_window/top5overall.csv")
length =dim(grand)[1]
results <- matrix(,length*10,7)
colnames(results) <-c("chr", "region","start","end","rsID","position","PIP")

for (i in 1:length) {
  
  chr = grand$chr[i]
  region = grand$region[i]
  
  mcmc <- read.table(paste0("/storage/nipm/kerimbae/sliding_window/chr",chr,"/c",region,".mcmc.txt"),header = T)
  
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

write.csv(results,"/storage/nipm/kerimbae/sliding_window/table2overallwithPIP.csv")

