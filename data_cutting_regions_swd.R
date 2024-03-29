# Benazir Rowe 
# 4/4/2018
# Cutting SSCCS chunks in mean genotype format
# accounting for the last leftover chunk

# bash command to start the R session on the remote cluster for high volume computations

# splitting into chunks recoded imputed data

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
    
    #five hundreds 500-1500,1500-2500, 2500-3500 etc. 
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
  
  #last region is composed by taking last 1000 SNPs of the chromosome
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







