# Benazir Rowe
# Summer 2019
# Create manhattan plot for MGS dataset

#input the number of regions per chromosome in SSCCS dataset

size <- matrix(,22,2)
size[,1]<- c(1:22)
size[,2]<- c(105,102,84,71,74,85,67,66,59,70,65,63,48,41,39,41,36,38,26,27,18,19)
size <- data.frame(size)
colnames(size)<- c("chr","chrlength")

############################ combining data by chrmosomes 1,3,5 chunks etc
for (chr in 1:22){
  chr = size$chr[chr]
  chrlength = size$chrlength[chr]
  
  #read first chunk
  placehold <- read.table(paste0("/storage/nipm/kerimbae/pimass/output/ssccs/chr",chr,"/c1.mcmc.txt"),header = TRUE)
  
  for (region in seq(3, chrlength, by = 2)){
    
    table <- read.table(paste0("/storage/nipm/kerimbae/pimass/output/ssccs/chr",chr,"/c",region,".mcmc.txt"),header = TRUE)
    placehold <- rbind(placehold,table)
    
  }
  
  if((chrlength %% 2) == 0) {
    print(paste(chrlength,"is Even"))
    
    
    table1 <- read.table(paste0("/storage/nipm/kerimbae/pimass/output/ssccs/chr",chr,"/c",chrlength,".mcmc.txt"),header = TRUE)
    table <- table1[501:1000,]
    placehold <- rbind(placehold,table)
    
  }
  
  dim(placehold)
  head(placehold)
  
  colnames(placehold)<-c("SNP","CHR","BP","P","beta")
  write.table(placehold,paste0("/storage/nipm/kerimbae/pimass/output/ssccs/chr",chr,"odd"))
}



### grand combine 22 chrs
grand <-  read.table(paste0("/storage/nipm/kerimbae/pimass/output/ssccs/chr1odd"),header=TRUE)

for (chr in 2:22){
  
  table <- read.table(paste0("/storage/nipm/kerimbae/pimass/output/ssccs/chr",chr,"odd"),header=TRUE)
  grand <- rbind(grand,table)
  dim(grand)
  print(tail(grand))
  
}

#combined dataset
write.table(grand,paste0("/storage/nipm/kerimbae/pimass/output/ssccs/grandodd"))

grand =read.table("/storage/nipm/kerimbae/pimass/output/ssccs/grandodd")

grand$beta <- 1-grand$P
colnames(grand)<-c("SNP","CHR","BP","Actual_P","P")

#combined dataset with 1-PIP instead of PIP for P
write.table(grand,paste0("/storage/nipm/kerimbae/pimass/output/ssccs/inversePnoprior"))

###PLOT 1.1 ############################################ basic Manhattan -log10 of (PIP) ###################################################

grand <-  read.table(paste0("/storage/nipm/kerimbae/pimass/output/ssccs/inversePnoprior"),header=TRUE)
range(grand$P)
-log10(0.10804)=0.9664154
-log10(0.99991)=3.908826e-05
quantile(q, c( .95,0.99))

library(qqman)
pdf(paste0("/storage/nipm/kerimbae/pimass/output/ssccs/figure1SWD.pdf"), width=10, height=5)

manhattan(grand,  ylim = c(0, 5), cex = 0.6, ylab=expression('log'[10]*'(1-PIP)'),
          cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F, 
          chrlabs = c(1:20, "21", "22"), 
          annotatePval = 0.5011872, 
          annotateTop = FALSE)

dev.off()

#scp kerimbae@cherry-creek.nscee.edu:/storage/nipm/kerimbae/pimass/output/ssccs/figure1SWD.pdf Desktop/

# Manhattan plot in .esp format

library(qqman)
grand <-  read.table(paste0("/storage/nipm/kerimbae/pimass/output/ssccs/inversePnoprior"),header=TRUE)

setEPS()

postscript("/storage/nipm/kerimbae/pimass/output/ssccs/figure4_7by5.eps",,width = 7,height = 5)

manhattan(grand,  ylim = c(0, 1), cex = 0.6, ylab=expression('log'[10]*'(1-PIP)'),
          cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F, 
          chrlabs = c(1:20, "21", "22"), 
          # annotatePval = 0.5 , 
          annotateTop = FALSE)

dev.off()


scp kerimbae@cherry-creek.nscee.edu:/storage/nipm/kerimbae/pimass/output/ssccs/figure4_7by5.eps Desktop/final_graphs




#Manhattan plot analysis
  
#subsetting required SNPs
newdata <- grand[order(-grand$P),]
q = newdata[1:20,]  
write.csv(q, "/storage/nipm/kerimbae/pimass/output/ssccs/top20SNPSSCCS.csv")

data <- read.csv("/storage/nipm/kerimbae/pimass/output/ssccs/summarySSCCS.csv")

# find the regions to which top 20 sNPs belong
for (i in 1:20){
  
  position = top20SNPsManh$BP[i]
  chr = top20SNPsManh$CHR[i]
  newdata <- data[which(position > data$start & position < data$end & data$chr == chr),]
  print(newdata)
  
}
  
