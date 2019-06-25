qsub -I -l ncpus=1,mem=2gb,cput=5:0:0 -l walltime=5:0:0 /bin/bash
module load intel intelmpi R
R

# retrieve positions and solve the mystery of chr 16

chr=16

#read position file
position.file = read.table(paste0("/storage/nipm/kerimbae/pimass/input/mgsInput/chr",chr,"posmgs.txt"))
position.file=data.frame(position.file)
dim(position.file)

#20170 by 3, which means we will have about 30 chunks

line = which(position.file[,1] == 'rs1386058')
position.file[line,]



