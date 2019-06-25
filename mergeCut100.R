qsub -I -l ncpus=1,mem=7gb,cput=5:0:0 -l walltime=5:0:0 /bin/bash
module load intel intelmpi R 
R

### data about dataset cut by 100 SNPs each no overlap

size <- matrix(,22,2)
size[,1]<- c(1:22)
size[,2]<- c(516,540,448,405,420,421,353,361,308,357,328,316,249,208,192,404,153,193,89,170,94,77)
size <- data.frame(size)
colnames(size)<-c("chr","chrlength")

## gather all regions by chromosomes

for (i in 1:22){
  chr=size$chr[i]
  chrlength=size$chrlength[i]
  
  #read first
  placehold <- read.table(paste0("/storage/nipm/kerimbae/pimass/output/newPprior/chr",chr,"/c1.mcmc.txt"),header = TRUE)
  
  for (i in 2:chrlength){
    
    table <- read.table(paste0("/storage/nipm/kerimbae/pimass/output/newPprior/chr",chr,"/c",i,".mcmc.txt"),header = TRUE)
    placehold <- rbind(placehold,table)
    
  }
  
  colnames(placehold)<-c("SNP","CHR","BP","P","beta")
  dim(placehold)
  head(placehold)
  
  
  placehold$beta <- 1-placehold$P
  colnames(placehold)<-c("SNP","CHR","BP","Actual_P","P")
  write.table(placehold,paste0("/storage/nipm/kerimbae/pimass/output/newPprior/chr",chr,"bySNP"))

}

####### combine 22

grand <-  read.table(paste0("/storage/nipm/kerimbae/pimass/output/newPprior/chr1bySNP"),header=TRUE)

for (chr in 2:22){
table <- read.table(paste0("/storage/nipm/kerimbae/pimass/output/newPprior/chr",chr,"bySNP"),header=TRUE)
grand <- rbind(grand,table)
dim(grand)
tail(grand)
}

write.table(grand,paste0("/storage/nipm/kerimbae/pimass/output/newPprior/grandBySNP"))

################################################# make a plot classic -log10

grand <-  read.table(paste0("/storage/nipm/kerimbae/pimass/output/newPprior/grandBySNP"),header=TRUE)

library(qqman)
pdf(paste0("/storage/nipm/kerimbae/pimass/output/newPprior/grandBySNPnewPprior.pdf"), width=10, height=5)

manhattan(grand, main = "Manhattan Plot", ylim = c(0, 4), cex = 0.6, 
          cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F, 
          chrlabs = c(1:20, "21", "22"))

dev.off()


#################### make a plot for 1.3

grand <-  read.table(paste0("/storage/nipm/kerimbae/pimass/output/newPprior/grandBySNP"),header=TRUE)

grand$P <- 10^(log(grand$P, 1.3))
write.table(grand,paste0("/storage/nipm/kerimbae/pimass/output/newPprior/grandBySNPlog13ready"))


library(qqman)
pdf(paste0("/storage/nipm/kerimbae/pimass/output/newPprior/grandBySNPlog13.pdf"), width=10, height=5)

manhattan(grand, main = "Manhattan Plot", ylim = c(0, 10), cex = 0.6, 
          cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F, 
          chrlabs = c(1:20, "21", "22"))

dev.off()



############################# plot by region



#read start positions
size <- matrix(,22,2)
size[,1]<- c(1:22)
size[,2]<- c(516,540,448,405,420,421,353,361,308,357,328,316,249,208,192,404,153,193,89,170,94,77)
size <- data.frame(size)
colnames(size)<-c("chr","chrlength")

#read start positions


for (i in 1:22){
    chr= size$chr[i]
    chrlength=size$chrlength[i]
    region=1
    position=1
    placehold <- read.table(paste0("/storage/nipm/kerimbae/pimass/output/newPprior/chr",chr,"/c1.mcmc.txt"),header = TRUE)
    placehold$beta <- 1-placehold$postc
    summa =sum(-log(placehold$beta))
    
    
    regions <- matrix(,chrlength,4)
    colnames(placehold)<-c("region","CHR","BP","P")
    regions[1,]<- c(region,chr,position,summa)
    head(regions)
    
      for (i in 2:chrlength){
        position=position+1
        placehold <- read.table(paste0("/storage/nipm/kerimbae/pimass/output/newPprior/chr",chr,"/c",i,".mcmc.txt"),header = TRUE)
        placehold$beta <- 1-placehold$postc
        summa =sum(-log(placehold$beta))
        regions[i,]<- c(i,chr,position,summa)
        
      }
    
    dim(regions)
    
    colnames(regions)<-c("SNP","CHR","BP","P")
    write.table(regions,paste0("/storage/nipm/kerimbae/pimass/output/newPprior/chrLog",chr,"byRegion"))
    tail(regions)

}

# combine 22 chrs in one data frame

grand <-  read.table(paste0("/storage/nipm/kerimbae/pimass/output/newPprior/chrLog1byRegion"),header=TRUE)

for (chr in 2:22){
table <- read.table(paste0("/storage/nipm/kerimbae/pimass/output/newPprior/chrLog",chr,"byRegion"),header=TRUE)
grand <- rbind(grand,table)
dim(grand)
tail(grand)

}
write.table(grand,paste0("/storage/nipm/kerimbae/pimass/output/newPprior/logByRegions"))


data <- read.table("logByRegions")
#pdf(paste0("LogByRegions100newPprior.pdf"), width=10, height=5)

plot(data$P  , pch=20,  col="blue",ylab = "sum(-log(1-PIP))",xlab="regionID", main= "Sum(-log(1-PIP) by region")
axis(1,at=data$CHR,pos=0,las=2)
#dev.off()  


######################## take a different log

grand$P <- 10^(log(grand$P, 1.3))



########################make plain sum by region

size <- matrix(,22,2)
size[,1]<- c(1:22)
size[,2]<- c(516,540,448,405,420,421,353,361,308,357,328,316,249,208,192,404,153,193,89,170,94,77)
size <- data.frame(size)
colnames(size)<-c("chr","chrlength")

#read start positions


for (i in 1:22){
  
  chr= size$chr[i]
  chrlength=size$chrlength[i]
  region=1
  position=1
  placehold <- read.table(paste0("/storage/nipm/kerimbae/pimass/output/newPprior/chr",chr,"/c1.mcmc.txt"),header = TRUE)
  summa =sum(placehold$postc)
  
  
  regions <- matrix(,chrlength,4)
  regions[1,]<- c(region,chr,position,summa)
  head(regions)
  
    for (i in 2:chrlength){
      position = position+1
      placehold <- read.table(paste0("/storage/nipm/kerimbae/pimass/output/newPprior/chr",chr,"/c",i,".mcmc.txt"),header = TRUE)
      summa <- sum(placehold$postc)
      regions[i,]<- c(i,chr,position,summa)
      
    }
  
  dim(regions)
  colnames(regions)<-c("SNP","CHR","BP","P")
  write.table(regions,paste0("/storage/nipm/kerimbae/pimass/output/newPprior/chr",chr,"SumbyRegion"))
  tail(regions)

}


# combine 22 chrs in one data frame

grand <-  read.table(paste0("/storage/nipm/kerimbae/pimass/output/newPprior/chr1SumbyRegion"),header=TRUE)

for (chr in 2:22){
  table <- read.table(paste0("/storage/nipm/kerimbae/pimass/output/newPprior/chr",chr,"SumbyRegion"),header=TRUE)
  grand <- rbind(grand,table)
  dim(grand)
  tail(grand)
}

write.table(grand,paste0("/storage/nipm/kerimbae/pimass/output/newPprior/SumByRegions"))

#plot Sum(???log(1???PIP))
data <- read.table("C:/Users/Benazir/Desktop/newPprior/logByRegions")
#pdf(paste0("LogByRegions100newPprior.pdf"), width=10, height=5)
plot(data$P  , pch=20,  col="blue",ylab = "Sum(???log(1???PIP))",xlab="regionID", main= "Sum(???log(1???PIP)) by region ")

#dev.off()  


# plot simple sum
data <- read.table("C:/Users/Benazir/Desktop/newPprior/SumByRegions")
#pdf(paste0("SumByRegions100newPprior.pdf"), width=10, height=5)

plot(data$P  , pch=20,  col="blue",ylab = "sum(PIP)",xlab="regionID", main= "Sum by region no h prior (cut by 100) ")

#dev.off()  

