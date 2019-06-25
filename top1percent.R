qsub -I -l ncpus=1,mem=10gb,cput=5:0:0 -l walltime=5:0:0 /bin/bash
module load intel intelmpi R
R

# based on the true phenotype run routput/gschri.txt want to subset top 5% of the regions.
# we do it based on their rank and check if percentile is >=95
# percentiles for top based on 10K

### data about dataset cut by 100 SNPs each no overlap

size <- matrix(,22,2)
size[,1]<- c(1:22)
size[,2]<- c(516,540,448,405,420,421,353,361,308,357,328,316,249,208,192,404,153,193,89,170,94,77)
size <- data.frame(size)
colnames(size)<-c("chr","chrlength")

for (i in 1:22){
 
   chr=size$chr[i]
   chrlength=size$chrlength[i]
   results <- matrix(,chrlength,5)
   
  for (region in 1:chrlength){
    
    mcmc <- read.table(paste0("/storage/nipm/kerimbae/pimass/output/cut100/chr",chr,"/c",region,".mcmc.txt"),header = T)
    summa = sum(mcmc$postc)
    start = mcmc$pos[1]
    end = mcmc$pos[dim(mcmc)[2]]
    
    results[region,] = rbind(chr, region,summa,start,end )
    
  }

  colnames(results) <-c("chr", "region","summa","start","end")
  
  results <- data.frame(results) 
  percentile <- ecdf(results$summa)(results$summa)
  ranks <- rank(-results$summa)
  grand <- cbind(results,percentile,ranks)
  
  write.csv(grand, paste0("/storage/nipm/kerimbae/pimass/output/cut100/summarychr",chr,".csv"),col.names = T,row.names = T,sep = " ")

}

############# phase 2: read top 1 percent regions and find the best SNP, create an array  
# name,location, chr,region

top1percent <- read.csv("/storage/nipm/kerimbae/pimass/output/cut100/top1percent.csv")

length = dim(top1percent)[1]

snps <- matrix(,length,5)

for (i in 1:length){
  
  chr = top1percent$chr[i]
  region = top1percent$region[i]
  
  mcmc <- read.table(paste0("/storage/nipm/kerimbae/pimass/output/cut100/chr",chr,"/c",region,".mcmc.txt"),header = T)
  index = which.max(mcmc$postc)
  snps[i,]<- rbind(as.character(mcmc$rs[index]),mcmc$pos[index],chr,region,mcmc$postc[index])
}
  
colnames(snps) <-c("rsID", "position","chr","region","value")

write.csv(snps, paste0("/storage/nipm/kerimbae/pimass/output/cut100/snps.csv"))


#### check chr 16 old result

chr=16
chrlength=39

results <- matrix(,chrlength,5)

for (region in 1:chrlength){
  
  mcmc <- read.table(paste0("/storage/nipm/kerimbae/pimass/output/chrgold",chr,"/c",region,".mcmc.txt"),header = T)
  summa = sum(mcmc$postc)
  start = mcmc$pos[1]
  end = mcmc$pos[dim(mcmc)[2]]
  
  results[region,] = rbind(chr, region,summa,start,end )
  
}

colnames(results) <-c("chr", "region","summa","start","end")

results <- data.frame(results) 
percentile <- ecdf(results$summa)(results$summa)
ranks <- rank(-results$summa)
grand <- cbind(results,percentile,ranks)

write.csv(grand, paste0("/storage/nipm/kerimbae/pimass/output/chrgold",chr,"/summary16.csv"))

