# Benazir Rowe 
# Spring 2018 UNLV 
# breaking up the data into overlapping regions for each chromosome for MGS dataset

# bash command to start the R session on the remote cluster for high volume computations

qsub -I -l ncpus=1,mem=5gb,cput=5:0:0 -l walltime=5:0:0 /bin/bash
module load intel intelmpi R 
R

# iterate over 22 chromosomes

for (q in 1:22){

path = paste("progs/pimass/pimass/input/singles/", q, "chr.txt",sep ="")
conn <- file(path,open="r")
lines <- readLines(conn)
d <- length(lines)
print(d)

#sliding window implementation
k = 1

    for (i in seq(1000, d, 1000)) {
      
      end <- i
      start <- i - 1000 + 1
      name <- paste("c", k, sep ="")
      assign(name, lines[start:end])
      
      write.table(get(name), file = paste("progs/pimass/pimass/input/chr",q,"/c", k ,".txt",sep =""), sep =",",quote = F, row.names = F, col.names = F)
      k = k + 1
      print(length(get(name)))
      start <- end - 500 + 1
      end <- end + 500
      
      name <- paste("c",k, sep ="")
      assign(name, lines[start:end])
      
      write.table(get(name), file = paste("progs/pimass/pimass/input/chr",q,"/c", k ,".txt",sep =""), sep =",",quote = F, row.names = F, col.names = F)
      k = k + 1
      print(length(get(name)))
    }

close (conn)

}






