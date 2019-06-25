qsub -I -l ncpus=1,mem=10gb,cput=5:0:0 -l walltime=5:0:0 /bin/bash
module load intel intelmpi R
R

# based on the true phenotype run routput/gschri.txt want to subset top 5% of the regions.
# we do it based on their rank and check if percentile is >=95
# percentiles for top based on 10K


chr = 1

pvals <- read.table(paste("routput/gschr",chr,".txt",sep=""),header = T)

percentile <- ecdf(pvals$sum)(pvals$sum)

# new<-sort(percentile)
# new
ranks <- rank(-pvals$sum)
ranks

grand<- cbind(pvals$sum,percentile,ranks)
grand

which(grand==1,arr.ind=TRUE)
which(grand==2,arr.ind=TRUE)
which(grand==3,arr.ind=TRUE)

which(grand==4,arr.ind=TRUE)


which(grand==5,arr.ind=TRUE)


which(grand==6,arr.ind=TRUE)

#### subset top 1% SNPs of the selected region 

table <- read.table("c6.mcmc.txt",header=T)
table.sorted = sort(table$postc, decreasing = T)

ranks<-rank(table$postc,ties.method = c("first"))
ranks

which(ranks==1,arr.ind=TRUE)
