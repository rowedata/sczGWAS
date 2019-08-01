# Benazir Rowe 
# Summer 2018 UNLV 
# create barplot of all regions and all chromosomes combined
# creat Manhattan plot

######################################################################################
#input the length of each chromosome MGS dataset

size <- matrix(,22,2)
size[,1] <- c(1:22)
size[,2] <- c(105,102,84,71,74,85,67,66,59,70,65,63,48,41,39,41,36,38,26,27,18,19)
size <- data.frame(size)
colnames(size)<- c("chr","chrlength")

#summarize each chromosomes mcmc results

for (k in 1:22){ #looping over chromosomes
 
  chr = size$chr[k]
  length = size$chrlength[k]  #chromosome length
  gs <- matrix(, length,7)    #matrix frame for chromosome's regions
  
   for (i in 1:length){ #looping over regions of chromosome
    
    infile <- paste("/storage/nipm/kerimbae/pimass/output/ssccs/chr",chr,"/c",i,".mcmc.txt", sep = "")
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
  write.table(gs, file =paste("/storage/nipm/kerimbae/pimass/output/ssccs/summarychr",chr,".txt",sep=""), row.names = F, col.names = T)
  
}


# combine 22 chrs results in one data frame

grand <-  read.table(paste0("/storage/nipm/kerimbae/pimass/output/ssccs/summarychr1.txt"),header=TRUE)

for (chr in 2:22){
  table <- read.table(paste0("/storage/nipm/kerimbae/pimass/output/ssccs/summarychr",chr,".txt"),header=TRUE)
  grand <- rbind(grand,table)
  print(dim(grand))
  print(tail(grand))
  
}
write.table(grand,paste0("/storage/nipm/kerimbae/pimass/output/ssccs/summaryswd.csv"),sep = ",",row.names = FALSE,
            col.names = TRUE) #save the combined file


# figure 1 barplot all regions based on sum PIP
allregions <- read.csv("summaryswd.csv",header=TRUE)

toString(allregions$CHR)
colvect <- matrix(,1244,1)

#indicator creation for alternating colors
for (i in 1:1244){
  if (allregions$chr[i] %%2 == 0){
    colvect[i]=1
  } else {
    colvect[i]=2
  }
}


allregions <- cbind(allregions,colvect) # add indicator vector to data frame 
library(ggplot2)

p <- ggplot(allregions, aes(x=chr, y=sum, colour=colvect))+
     geom_bar(stat = "identity", position = "dodge2",show.legend = FALSE)
  p + scale_fill_manual(values = c("blue4", "darkorange"))
  
  p + scale_x_continuous(breaks=seq(1,22,1))+labs(y="SUM(PIP)",x="Chromosome")
 



#manhattan plot
  
  grand <-  read.table(paste0("/storage/nipm/kerimbae/pimass/output/ssccs/grandBySNP"),header=TRUE)
  
  library(qqman)
  pdf(paste0("/storage/nipm/kerimbae/pimass/output/ssccs/manhattanSWD.pdf"), width=10, height=5)
  
  manhattan(grand, main = "Manhattan Plot", ylim = c(0, 4), cex = 0.6, 
            cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F, 
            chrlabs = c(1:20, "21", "22"),xlab="Chromosome",ylab=expression('log'[10]*'(1-PIP)'))
  
  dev.off()


#locally each chromosome separately

mcmc <- read.table(paste0("C:/Users/Benazir/Desktop/ssccs/resultchr",chr,".txt"),header=T)
length=dim(mcmc)[1]
barplot(mcmc$sum, main=paste("Chromosome ", chr," Sum of Posterior Inclusion Probabilities"), xlab="Group number", names.arg= c(1:length))
abline(h = quantile(as.numeric(mcmc$sum), c(.9), na.rm=T), lty=2)
abline(h = quantile(as.numeric(mcmc$sum), c(.75), na.rm=T), col=c("red"))
legend(1, 1, c("90th percentile", "75th percentile"), col = c("black", "red"),lty = c(2, 1),lwd=c(1,1),bg="gray95")

# stat analysis overlap mgs ssccs

newdata <- grand[order(-sum),]
top5regions <- newdata[1:63,]

#first
i=1
chr = top5regions$chr[i]
region =top5regions$region[i]

file <- read.table(paste0("/storage/nipm/kerimbae/pimass/output/ssccs/chr",chr,"/c",region,".mcmc.txt"),header=T)
orderedmcmc <- file[order(-file$postc),]
table <- orderedmcmc[1:10, ]


vect <- matrix(region, 10, 1)
q <- as.character(table$rs)
newtable <- cbind(as.numeric(table$chr),vect,min(file$pos),max(file$pos),q,table$pos,table$postc)
supplTable <- newtable



for (i in 2:63) {
  
  chr = top5regions$chr[i]
  region =top5regions$region[i]

  file <- read.table(paste0("/storage/nipm/kerimbae/pimass/output/ssccs/chr",chr,"/c",region,".mcmc.txt"),header=T)
  orderedmcmc <- file[order(-file$postc),]
  table <- orderedmcmc[1:10,]
  vect <- matrix(region,10,1)
  q <-as.character(table$rs)
  
  newtable <- cbind(as.numeric(table$chr),vect,min(file$pos),max(file$pos),q,table$pos,table$postc)
  supplTable <- rbind(supplTable, newtable)
  dim(supplTable)
  tail(supplTable)
  
}

colnames(supplTable)<-c("chromosome","region","start", "end","rsID","location","PIP")



write.csv(supplTable,"/storage/nipm/kerimbae/pimass/output/ssccs/supplTableSWD.csv", quote=FALSE,row.names = FALSE)






