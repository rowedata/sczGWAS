qsub -I -l ncpus=1,mem=10gb,cput=5:0:0 -l walltime=5:0:0 /bin/bash
module load intel intelmpi R
R

nineK.dir <- "/storage/nipm/kerimbae/pvalue/"
oneK.dir <- ""
hundredK.dir <- "/storage/nipm/kerimbae/pimass/output/"

######################################################################################
#change input ONLY here
# ALL RUNS in CHERRY CREEK, FILES ARE IMPORTED AFTER PROCESSING

#########################################

chr = 4
chunk= 74
output <- read.table(paste("routput/gschr",chr,".txt",sep=""),header = T)
