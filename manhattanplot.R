### RUN 1################################################

qsub -I -l ncpus=1,mem=7gb,cput=5:0:0 -l walltime=5:0:0 /bin/bash
module load intel intelmpi R 
R


# plot type 1: merging odds
# If odd number of chunks great just merge. If even cut and add a half of the last one

# hand input the chromosome size

size <- matrix(,22,2)
size[,1]<- c(1:22)
size[,2]<- c(103,107,89,80,83,84,70,72,61,71,65,63,49,41,38,39,30,38,17,33,18,15)
size <- data.frame(size)
colnames(size)<-c("chr","chrlength")

############################ combining data by chrmosomes 1,3,5 chunks etc
for (chr in 1:22){
  chr = size$chr[chr]
  chrlength = size$chrlength[chr]
  
  #read first chunk
  placehold <- read.table(paste0("/storage/nipm/kerimbae/sliding_window/chr",chr,"/c1.mcmc.txt"),header = TRUE)
  
  for (region in seq(3, chrlength, by = 2)){
    
    table <- read.table(paste0("/storage/nipm/kerimbae/sliding_window/chr",chr,"/c",region,".mcmc.txt"),header = TRUE)
    placehold <- rbind(placehold,table)
    
  }
  
  if((chrlength %% 2) == 0) {
    print(paste(chrlength,"is Even"))
    
    
    table1 <- read.table(paste0("/storage/nipm/kerimbae/sliding_window/chr",chr,"/c",chrlength,".mcmc.txt"),header = TRUE)
    table <- table1[501:1000,]
    placehold <- rbind(placehold,table)
    
  }
  
  dim(placehold)
  head(placehold)
  
  colnames(placehold)<-c("SNP","CHR","BP","P","beta")
  write.table(placehold,paste0("/storage/nipm/kerimbae/sliding_window/chr",chr,"odd"))
}



### grand combine 22 chrs
grand <-  read.table(paste0("/storage/nipm/kerimbae/sliding_window/chr1odd"),header=TRUE)

for (chr in 2:22){
  
  table <- read.table(paste0("/storage/nipm/kerimbae/sliding_window/chr",chr,"odd"),header=TRUE)
  grand <- rbind(grand,table)
  dim(grand)
  print(tail(grand))
  
}

#combined dataset
write.table(grand,paste0("/storage/nipm/kerimbae/sliding_window/grandodd"))

grand =read.table("/storage/nipm/kerimbae/sliding_window/grandodd")

grand$beta <- 1-grand$P
colnames(grand)<-c("SNP","CHR","BP","Actual_P","P")

#combined dataset with 1-PIP instead of PIP for P
write.table(grand,paste0("/storage/nipm/kerimbae/sliding_window/inversePnoprior"))

###PLOT 1.1 ############################################ basic Manhattan -log10 of (PIP) ###################################################


qsub -I -l ncpus=1,mem=7gb,cput=5:0:0 -l walltime=5:0:0 /bin/bash
module load intel intelmpi R 
R

grand <-  read.table(paste0("/storage/nipm/kerimbae/sliding_window/inversePnoprior"),header=TRUE)
range(grand$P)
-log10(0.10804)=0.9664154
-log10(0.99991)=3.908826e-05
quantile(q, c( .95,0.99))

library(qqman)
pdf(paste0("/storage/nipm/kerimbae/sliding_window/figure1_4.pdf"), width=10, height=5)

manhattan(grand,  ylim = c(0, 1), cex = 0.6, ylab=expression('log'[10]*'(1-PIP)'),
          cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F, 
          chrlabs = c(1:20, "21", "22"), annotatePval = 0.5011872, annotateTop = FALSE)

dev.off()


scp kerimbae@cherry-creek.nscee.edu:/storage/nipm/kerimbae/sliding_window/figure1_4.pdf Desktop/

# PLOT 1.2 ###############################################################################################################
#create new dataset, replace P for new value
grand <-  read.table(paste0("/storage/nipm/kerimbae/sliding_window/inversePnoprior"),header=TRUE)
newgrandP = 10^(log(grand$P, base=1.3))
grand$P= newgrandP
head(grand)

write.table(grand,paste0("/storage/nipm/kerimbae/sliding_window/inversePnopriorlog13"))

##########
grand <-  read.table(paste0("/storage/nipm/kerimbae/sliding_window/inversePnopriorlog13"),header=TRUE)
library(qqman)
pdf(paste0("/storage/nipm/kerimbae/sliding_window/log13noprior.pdf"), width=10, height=5)

manhattan(grand, main = "Manhattan Plot", ylim = c(0, 30), cex = 0.6, 
          cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F, 
          chrlabs = c(1:20, "21", "22"))

dev.off()

#download the plot from cherry-creek
scp kerimbae@cherry-creek.nscee.edu:/storage/nipm/kerimbae/sliding_window/log13noprior.pdf Desktop/slidingW



######PLOT 1.3###############
############################################################################

for (chr in 1:22){
  #set start values for chromosomes
  chr = size$chr[chr]
  chrlength = size$chrlength[chr]
  
  placehold <- read.table(paste0("/storage/nipm/kerimbae/sliding_window/chr",chr,"/c1.mcmc.txt"),header = TRUE)
  
  q =  1-placehold$postc
  logsum = sum(-log10(q))
  summa =sum(placehold$postc)
  
  
  regions <- matrix(,chrlength,4)
  colnames(regions)<-c("CHR","Region","SUM_PIP","LOGSUM")
  regions[1,]<- c(chr ,1, summa, logsum)
  head(regions)
  
  for (region in 2:chrlength){
    
    placehold <- read.table(paste0("/storage/nipm/kerimbae/sliding_window/chr",chr,"/c",region,".mcmc.txt"),header = TRUE)
    q =  1-placehold$postc
    logsum = sum(-log10(q))
    summa =sum(placehold$postc)
    regions[region,]<- c(chr,region,summa,logsum)
    
  }
  
  write.table(regions,paste0("/storage/nipm/kerimbae/sliding_window/chr",chr,"regions"))
  print(dim(regions)) 
  print (tail(regions))
  
}

# combine 22 chrs in one data frame

grand <-  read.table(paste0("/storage/nipm/kerimbae/sliding_window/chr1regions"),header=TRUE)

for (chr in 2:22){
  table <- read.table(paste0("/storage/nipm/kerimbae/sliding_window/chr",chr,"regions"),header=TRUE)
  grand <- rbind(grand,table)
  dim(grand)
  tail(grand)
}
write.table(grand,paste0("/storage/nipm/kerimbae/sliding_window/allregions"))



######################################################################################## Log plot
allregions <-  read.table(paste0("/storage/nipm/kerimbae/sliding_window/allregions"),header=TRUE)

pdf(paste0("/storage/nipm/kerimbae/sliding_window/figure2.pdf"), width=10, height=5)

barplot(allregions$SUM_PIP  , pch=20,  col="blue",ylab = "SUM(PIP)", xlab="region", main= "Sum by region  (cut by 1000) ")
#plot(allregions$LOGSUM  , pch=20,  col="blue",ylab = "SUM(-log10(1-PIP))", xlab="regionID", main= "SUM(-log10(1-PIP)) by region (cut by 1000) ")

dev.off()

scp kerimbae@cherry-creek.nscee.edu:/storage/nipm/kerimbae/sliding_window/figure2.pdf Desktop/
# pdf are in folder slidingW

# on local R to play with graph
allregions <-  read.table("allregions",header=TRUE)
barplot(allregions$SUM_PIP  , pch=20,  col="blue",ylab = "SUM(PIP)", xlab="region", main= "Sum by region  (cut by 1000) ")

library(ggplot2)
p  <- ggplot(allregions,aes(x=CHR,y=SUM_PIP),color=CHR) +
  geom_bar(stat = "identity", aes(fill = SUM_PIP), position = "dodge2")+
  scale_fill_manual(breaks = c("2", "1", "0.5"), 
                    values=c("red", "blue", "red", "blue","red", "blue","red", "blue",
                             "red", "blue","red", "blue","red", "blue","red", "blue",
                             "red", "blue","red", "blue","red", "blue","red", "blue",
                             "red", "blue","red", "blue"))
  
p






####END of PLOTS#####
###################making decision whether sliding window is worth it

chr=16

#read first
c1 <- read.table(paste0("/storage/nipm/kerimbae/sliding_window/chr",chr,"/c1.mcmc.txt"),header = TRUE)
c2 <- read.table(paste0("/storage/nipm/kerimbae/sliding_window/chr",chr,"/c2.mcmc.txt"),header = TRUE)

q =c1$postc[501:1000]
p= c2$postc[1:500]
diff = q-p
summary(diff)
summary(p)
summary(q)

pdf(paste0("/storage/nipm/kerimbae/sliding_window/diff161.pdf"), width=10, height=5)
par(mfrow=c(3,1))
barplot(p, main="0.000610 0.001020 0.001285 0.002961 0.001830 0.062440")
barplot(q,main="0.0004700 0.0008775 0.0011500 0.0025850 0.0016700 0.0461200")
barplot(diff, main="-0.0163200 -0.0003800 -0.0001500 -0.0003759  0.0001000  0.0057200")

dev.off()


###################################################################
### Try the result for 100SNPs ###
## take the grand SNP file and cut by position.

qsub -I -l ncpus=1,mem=7gb,cput=5:0:0 -l walltime=5:0:0 /bin/bash
module load intel intelmpi R 
R

size <- matrix(,22,2)
size[,1]<- c(1:22)
size[,2]<- c(103,107,89,80,83,84,70,72,61,71,65,63,49,41,38,39,30,38,17,33,18,15)
size <- data.frame(size)
colnames(size)<-c("chr","chrlength")

############################# sum PIP by 100
grand <-  read.table(paste0("/storage/nipm/kerimbae/sliding_window/grandodd"),header=TRUE)



length = dim(grand)[1]

q <- matrix(,6440,1)

k = 1
for (i in seq(100, length, by = 100)){
  
  start=i-100+1
  end=i
  summ= sum(grand$P[start:end])
  q[k]=summ
  k = k + 1
}


write.table(q,paste0("/storage/nipm/kerimbae/sliding_window/1000by100"))
data<- read.table("1000by100")
plot(data$V1, pch=20,  col="blue",ylab = "Sum(PIP)",xlab="new regionID", main= "Sum(PIP) for sliding window run cut by 100")

##########Sum(-log10(1-PIP))for cut by 100
grand <-  read.table(paste0("/storage/nipm/kerimbae/sliding_window/grandodd"),header=TRUE)

length = dim(grand)[1]
q <- matrix(,6440,1)
c <-matrix(1,100,1)
k = 1

for (i in seq(100, length, by = 100)){
  
  start=i-100+1
  end=i
  summ= -log10(sum(c-grand$P[start:end]))
  q[k]=summ
  k = k + 1
}


write.table(q,paste0("/storage/nipm/kerimbae/sliding_window/log1000by100"))
data<- read.table("log1000by100")
plot(data$V1, pch=20,  col="blue",ylab = "Sum(-log10(1-PIP))",xlab="new regionID", main= "Sum(-log10(1-PIP)) for sliding window run cut by 100")






#########################################33


for (i in 1:22){
  chr=size$chr[i]
  chrlength=size$chrlength[i]
  
  region=1
  position=1
  placehold <- read.table(paste0("/storage/nipm/kerimbae/sliding_window/chr",chr,"/c1.mcmc.txt"),header = TRUE)
  summa =sum(placehold$postc)
  
  
  regions <- matrix(,chrlength,4)
  colnames(placehold)<-c("region","CHR","BP","P")
  regions[1,]<- c(region,chr,position,summa)
  head(regions)
  
  for (i in 2:chrlength){
    position = position+1
    table <- read.table(paste0("/storage/nipm/kerimbae/sliding_window/chr",chr,"/c",i,".mcmc.txt"),header = TRUE)
    summa <- sum(table$postc)
    regions[i,]<- c(i,chr,position,summa)
    
  }
  
  dim(regions)
  colnames(regions)<-c("SNP","CHR","BP","P")
  write.table(regions,paste0("/storage/nipm/kerimbae/sliding_window/sumchr",chr,"regions"))
  print(tail(regions))
}





# combine 22 chrs in one data frame

grand <-  read.table(paste0("/storage/nipm/kerimbae/sliding_window/sumchr1regions"),header=TRUE)

for (chr in 2:22){
  table <- read.table(paste0("/storage/nipm/kerimbae/sliding_window/sumchr",chr,"regions"),header=TRUE)
  grand <- rbind(grand,table)
  dim(grand)
  tail(grand)
  
}

write.table(grand, "/storage/nipm/kerimbae/sliding_window/dataRegions")
data<- read.table("dataRegions")
plot(data$P, pch=20,  col="blue",ylab = "Sum(PIP)",xlab="regionID", main= "Sum(PIP) for sliding window")

###########calculate size of the region





#sum PIP for 1000
#sum PIP for 1000

for (i in 1:22){
  chr=size$chr[i]
  chrlength=size$chrlength[i]
  
  region=1
  position=1
  placehold <- read.table(paste0("/storage/nipm/kerimbae/sliding_window/chr",chr,"/c1.mcmc.txt"),header = TRUE)
  summa =sum(placehold$postc)
  range=placehold$pos[1000]-placehold$pos[1]
  
  regions <- matrix(,chrlength,5)
  regions[1,]<- c(region,chr,position,summa,range)
  head(regions)
  
  for (i in 2:chrlength){
    position = position+1
    table <- read.table(paste0("/storage/nipm/kerimbae/sliding_window/chr",chr,"/c",i,".mcmc.txt"),header = TRUE)
    summa <- sum(table$postc)
    range=table$pos[1000]-table$pos[1]
    regions[i,]<- c(i,chr,position,summa,range)
    
  }
  
  dim(regions)
  colnames(regions)<-c("SNP","CHR","BP","P","Range")
  write.table(regions,paste0("/storage/nipm/kerimbae/sliding_window/sumchr",chr,"Range"))
  print(tail(regions))
}


grand <-  read.table(paste0("/storage/nipm/kerimbae/sliding_window/sumchr1Range"),header=TRUE)

for (chr in 2:22){
  table <- read.table(paste0("/storage/nipm/kerimbae/sliding_window/sumchr",chr,"Range"),header=TRUE)
  grand <- rbind(grand,table)
  dim(grand)
  tail(grand)
  
}

write.table(grand, "/storage/nipm/kerimbae/sliding_window/dataRange")












