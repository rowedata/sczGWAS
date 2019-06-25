# Benazir Rowe 
# Recoding swedish dataset from bimbam 1,2 to chunks in mean genotype format
#4/4/2018
# accounting for the last leftover chunk

#code to start R on cluster
qsub -I -l ncpus=1,mem=5gb,cput=5:0:0 -l walltime=5:0:0 /bin/bash
module load intel intelmpi R 
R

#actual code for splitting into chunks recoded imputed data

for (q in 1:22){
  
  #k is a counter for the number of chunks, start from 1 for each chromosome
  k = 1
  
  path = paste("swdrecode/chr", q, ".txt",sep ="")
  conn <- file(path,open="r")
  lines <- readLines(conn)
  d <- length(lines)
  print(d)
  
  #open connection for log file
  con <- file(paste("swdrecode/swdrecode", q,".log", sep ="")) #save at the source
  sink(con, append=TRUE, type="output")

  for (i in seq(1000, d, 1000)) {
    
    #zeros 
    end <- i
    start <- i - 1000 + 1
    name <- paste("c", k, sep ="")
    assign(name, lines[start:end])
    
    write.table(get(name), file = paste("progs/pimass/pimass/input/swd/chr",q,"/c", k ,".txt",sep =""), sep =",",quote = F, row.names = F, col.names = F)
    
    print(k)
    print(start)
    print(end)
    k = k + 1
    
    #five hundreds  
    end <- end + 500
    if (end < d){
      name <- paste("c",k, sep ="")
      assign(name, lines[start:end])
      
      write.table(get(name), file = paste("progs/pimass/pimass/input/swd/chr",q,"/c", k ,".txt",sep =""), sep =",",quote = F, row.names = F, col.names = F)
      
      
      print(k)
      print(start)
      print(end)
      k = k + 1
    }
  }
  
  #last chunk leftover
  name <- paste("c", k, sep ="")
  
  assign(name, lines[(d - 1000 + 1):d])
  write.table(get(name), file = paste("progs/pimass/pimass/input/swd/chr",q,"/c", k ,".txt",sep =""), sep =",",quote = F, row.names = F, col.names = F)
  print(k)
  print(d - 1000 + 1)
  print(d)
  
  close (conn)
  
  # close the output log file connenction
  sink()
  close (con)
}







