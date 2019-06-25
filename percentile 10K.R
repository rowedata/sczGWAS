qsub -I -l ncpus=1,mem=1gb,cput=5:0:0 -l walltime=5:0:0 /bin/bash
module load intel intelmpi R
R

# percentiles for top based on 10K

chr = 16
chunk = 16
pvals <- read.table(paste("routput/gschr",chr,".txt",sep=""),header = T)

ranks <- rank(-pvals$sum)
ranks[chunk]

percentile <- ecdf(pvals$sum)(pvals$sum)
percentile[chunk]
