# Benazir Rowe
# Fall 2018 UNLV
# simulating permuted phenotype
# Inputs: original phenotype vector, count the number of 1s
# Output: 10,000 files containing permuted phenotypes to use in permutation test

# read original phenotype
pheno <- read.table("progs/pimass/pimass/input/mgsInput/phenoscz.txt")
tp <- t(pheno)

#multiple time permutation
for (i in 90001:100000)
{

  vect <- sample (tp) #permutation procedure
  write.table(vect,paste("/storage/nipm/kerimbae/pimass/input/pheno100K/ph",i,".txt", sep=""), row.names = FALSE, col.names = FALSE )
  
}



