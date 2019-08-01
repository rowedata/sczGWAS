# Benazir Rowe 
# Summer 2019 UNLV 
# Statistical analysis of the overlap between MGS and SSCCS datasets
# Overlap is defined based on various distances



# bash command to start the R session on the remote cluster for high volume computations
qsub -I -l ncpus=1,mem=5gb,cput=5:0:0 -l walltime=5:0:0 /bin/bash
module load intel intelmpi R
R

# number of regions for each chromosome of the SSCCS dataset, manual input
size <- matrix(,22,2)
size[,1]<- c(1:22)
size[,2]<- c(103,107,89,80,83,84,70,72,61,71,65,63,49,41,38,39,30,38,17,33,18,15)
size <- data.frame(size)
colnames(size)<-c("chr","chrlength")


# stat analysis overlap mgs ssccs
grand <- read.csv("/storage/nipm/kerimbae/sliding_window/summaryMGS.csv")

newdata <- grand[order(-grand$sum),]

top5regions <- newdata[1:64,]

#first
i = 1
chr = top5regions$chr[i]
region =top5regions$region[i]

file <- read.table(paste0("/storage/nipm/kerimbae/sliding_window/chr",chr,"/c",region,".mcmc.txt"),header=T)
orderedmcmc <- file[order(-file$postc),]
table <- orderedmcmc[1:10,]


vect <- matrix(region,10,1)
q <-as.character(table$rs)
newtable <- cbind(as.numeric(table$chr),vect,min(file$pos),max(file$pos),q,table$pos,table$postc)
supplTable <- newtable



for (i in 2:64) {
  
  chr = top5regions$chr[i]
  region =top5regions$region[i]
  
  
  file <- read.table(paste0("/storage/nipm/kerimbae/sliding_window/chr",chr,"/c",region,".mcmc.txt"),header=T)
  orderedmcmc <- file[order(-file$postc),] #order by sum PIP
  table <- orderedmcmc[1:10,]              #subset top10
  vect <- matrix(region,10,1)
  q <-as.character(table$rs)
  
  newtable <- cbind(as.numeric(table$chr),vect,min(file$pos),max(file$pos),q,table$pos,table$postc)
  supplTable <- rbind(supplTable,newtable)
  dim(supplTable)
  tail(supplTable)
  
}

colnames(supplTable)<-c("chromosome","region","start", "end","rsID","location","PIP")

write.csv(supplTable,"/storage/nipm/kerimbae/sliding_window/supplTableMGS.csv", quote=FALSE,row.names = FALSE)


# statistical analysis of the overlap between MGS and SSCCS 
 
setwd("C:/Users/Benazir/Documents/")

MGS <- read.csv("PAPER ONE/supplTableMGS.csv")
SSCCS <- read.csv("PAPER ONE/supplTableSWD.csv")

#subset the same chromosomes

result <- result[FALSE,]

for (i in 1:22){
  chrMGS = subset(MGS, MGS$chromosome==i)
  chrSSCCS = subset(SSCCS, SSCCS$chromosome==i)
  
  lengthMGS = dim(chrMGS)[1]
  lengthSSCCS = dim(chrSSCCS)[1]
  if (lengthMGS==0 |lengthSSCCS==0){
    next
  }
  for (m in 1:lengthMGS){
    
    for (s in 1:lengthSSCCS){
      
      distance =abs(chrMGS$location[m]-chrSSCCS$location[s])
      if (distance < 100000){
        print(distance)
        
        newfile =cbind(chrMGS[m,],chrSSCCS[s,])
        result <- rbind(result,newfile)
      }
      
    }
  }
  
  
}

result
dim(result)[1]
length =dim(result)[1]
round((dim(result)[1]/640) ,2)


for (i in 1:(length-1)){
  if (result[i,5]==result[(i+1),5] & result[i,12]==result[(i+1),12]){
    print (result[i:(i+1),])
  }
  
}







