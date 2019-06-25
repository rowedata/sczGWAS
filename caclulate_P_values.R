qsub -I -l ncpus=1,mem=1gb,cput=5:0:0 -l walltime=5:0:0 /bin/bash
module load intel intelmpi R
R

#nineK.dir <- "/storage/nipm/kerimbae/pvalue/"
#oneK.dir <- ""
hundredK.dir <- "/storage/nipm/kerimbae/pimass/output/"

######################################################################################
#change input ONLY here
# ALL RUNS in CHERRY CREEK, FILES ARE IMPORTED AFTER PROCESSING

#########################################

chr = 14
chunk= 33
output <- read.table(paste("routput/gschr",chr,".txt",sep=""),header = T)



start=20001
end =30000

folder.num=3



#################################################################
# single pvalue calculation, less than 100 (by hand)

k = 0
for (i in start:end){
  
  
  infile <- paste(hundredK.dir,chr,chunk,"/tenK",folder.num,"/",chr,"chr",chunk,"run",i,".mcmc.txt", sep = "")
  mcmc <- read.table(infile, head = T, stringsAsFactors = FALSE)
  
  if( sum(mcmc$postc) > output[chunk,1])
  {
    k = k + 1
    print(i)
  }
  
}

k






