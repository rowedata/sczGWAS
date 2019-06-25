#Benazir Rowe: simulating permuted phenotype
# Inputs: original phenotype vector, count the number of 1s
# Output: 10,000 files containing permuted phenotypes to use in permutation test

qsub -I -l ncpus=1,mem=10gb,cput=5:0:0 -l walltime=5:0:0 /bin/bash
module load intel intelmpi R
R


#one time read
pheno <- read.table("progs/pimass/pimass/input/mgsInput/phenoscz.txt")
tp <- t(pheno)

#multiple time permutation


for (i in 90001:100000)
{

  vect <- sample (tp)
  write.table(vect,paste("/storage/nipm/kerimbae/pimass/input/pheno100K/ph",i,".txt", sep=""), row.names = FALSE, col.names = FALSE )
  
}



