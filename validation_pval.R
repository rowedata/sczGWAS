qsub -I -l ncpus=1,mem=2gb,cput=5:0:0 -l walltime=5:0:0 /bin/bash
module load intel intelmpi R
R


validation.dir <- "progs/pimass/pimass/output/validation"

######################################################################################
#running 10,000 validations on last 18 chunks, calculating results here
# progs/pimass/pimass/input/newpbs/validation_true.pbs - pbs file to run jobs



chr = 4
chunk=74

start=1
end = 1000
##########################
mcmc_gold <- read.table(paste0(validation.dir,"/true",chr,chunk,".mcmc.txt"),header = T)
gold = mean(mcmc_gold$postc)


# single pvalue calculation, less than 100 (by hand)

k = 0
for (i in start:end){
  
  infile <- paste0(validation.dir,"/chr",chr,chunk,"/c",chr,chunk,"run",i,".mcmc.txt")
  mcmc <- read.table(infile, head = T, stringsAsFactors = FALSE)
  
  if( mean(mcmc$postc) > gold)
  {
    k = k + 1
    print(i)
  }
  
}

k



#after moving to /storage

start=895

for (i in start:end){
  
  infile <- paste0("/storage/nipm/kerimbae/pimass/output/validation/chr",chr,chunk,"/c",chr,chunk,"run",i,".mcmc.txt")
  mcmc <- read.table(infile, head = T, stringsAsFactors = FALSE)
  
  if( mean(mcmc$postc) > gold)
  {
    k = k + 1
    print(i)
  }
  
}

k

