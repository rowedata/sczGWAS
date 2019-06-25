qsub -I -l ncpus=1,mem=5gb,cput=5:0:0 -l walltime=5:0:0 /bin/bash
module load intel intelmpi R 
R


# plot type 1: merging odds
# If odd number of chunks great just merge. If even cut and add a half of the last one

# hand input the chromosome size
winners <- matrix(,18,2)
winners[,1]<- c(4,8,9,9,9,14,14,15,18,19,20,1,3,5,13,15,7,16)
winners[,2]<- c(74,30,23,24,32,6,33,29,15,5,1,36,15,34,41,28,58,16)
winners <- data.frame(winners)
colnames(winners)<- c("chr","region")

###

size <- matrix(,22,2)
size[,1]<- c(1:22)
size[,2]<- c(103,107,89,80,83,84,70,72,61,71,65,63,49,41,38,39,30,38,17,33,18,15)
size <- data.frame(size)
colnames(size)<-c("chr","chrlength")




#summarize each chromosomes mcmc results

for (k in 1:22){ #looping over chromosomes
  
  chr = size$chr[k]
  length = size$chrlength[k]  #chromosome length
  gs <- matrix(, length,7)    #matrix frame for chromosome's regions
  
  for (i in 1:length){ #looping over regions of chromosome
    
    infile <- paste("/storage/nipm/kerimbae/sliding_window/chr",chr,"/c",i,".mcmc.txt", sep = "")
    mcmc <- read.table(infile, head = T, stringsAsFactors = FALSE)
    
    gs[i,1] <- sum(mcmc$postc)
    gs[i,2] <- mean(mcmc$postc)
    gs[i,3] <- max(mcmc$postc)
    gs[i,4] <- min(mcmc$pos)
    gs[i,5] <- max(mcmc$pos)
    gs[i,6] <- chr
    gs[i,7] <- i
    
  }
  
  colnames(gs) = c("sum","mean","max", "start","end","chr","region")
  write.table(gs, file =paste("/storage/nipm/kerimbae/sliding_window/summarychr",chr,".txt",sep=""), row.names = F, col.names = T)
  
}


# combine 22 chrs results in one data frame

grand <-  read.table(paste0("/storage/nipm/kerimbae/sliding_window/summarychr1.txt"),header=TRUE)

for (chr in 2:22){
  table <- read.table(paste0("/storage/nipm/kerimbae/sliding_window/summarychr",chr,".txt"),header=TRUE)
  grand <- rbind(grand,table)
  print(dim(grand))
  print(tail(grand))
  
}
write.table(grand,paste0("/storage/nipm/kerimbae/sliding_window/summaryMGS.csv"),sep = ",",row.names = FALSE,
            col.names = TRUE)

newdata <- grand[order(-sum),]

ranks <- matrix(1:1266,1266,1)
withranks <- cbind(newdata,ranks)


for (i in 1:18){
  rankdata <- withranks[ which(withranks$chr==winners$chr[i] & withranks$region==winners$region[i]), ]
  print(rankdata)
  
  write.table(rankdata,file="/storage/nipm/kerimbae/sliding_window/ranksMGSpermutation.txt",append=TRUE)
  
}



