# Benazir Rowe
# Summer 2019 UNLV
# Generation of manhattan plot (single SNP PIP) and barplot for regions (sums of 1000 SNPs)

##Figure 1 ############################################ basic Manhattan -log10 of (PIP) ############


grand <-  read.table("inversePnoprior", header=TRUE)

range(grand$P)
-log10(0.10804)=0.9664154
-log10(0.99991)=3.908826e-05
quantile(q, c( .95,0.99))

library(qqman)
pdf(paste0("figure2_Feb5.pdf"), width=10, height=5)

manhattan(grand,  ylim = c(0, 1), cex = 0.6, ylab=expression('log'[10]*'(1-PIP)'),
          cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F, 
          chrlabs = c(1:20, "21", "22"))

dev.off()


#scp kerimbae@cherry-creek.nscee.edu:/storage/nipm/kerimbae/sliding_window/figure1_4.pdf Desktop/

#### Figure 2 ############################################ Log plot


allregions <-  read.table(paste0("/storage/nipm/kerimbae/sliding_window/allregions"),header=TRUE)

pdf(paste0("/storage/nipm/kerimbae/sliding_window/figure2.pdf"), width=10, height=5)

barplot(allregions$SUM_PIP  , pch=20,  col="blue",ylab = "SUM(PIP)", xlab="region", main= "Sum by region  (cut by 1000) ")

dev.off()


#### CODE FOR Figure 2

allregions <-  read.table("allregions",header=TRUE)
#barplot(allregions$SUM_PIP  , pch=20,  col="blue",ylab = "SUM(PIP)", xlab="region", main= "Sum by region  (cut by 1000) ")
toString(allregions$CHR)

#create a vector to address the alternating color
colvect <- matrix(,1266,1)

#new
for (i in 1:1266){
  if (allregions$CHR[i] %%2==0){
    colvect[i]=1
  } else {
    colvect[i]=2
  }
}


allregions <- cbind(allregions,colvect)
library(ggplot2)

p  <- ggplot(allregions, aes(x=CHR,y=SUM_PIP, fill=colvect))+ 
  geom_bar(stat = "identity",  position = "dodge2",show.legend = FALSE)
  
p + scale_x_continuous(breaks=seq(1,22,1))+labs(y="SUM(PIP)",x="Chromosome")






