# Benazir Rowe 
# Summer 2018 UNLV 
# Cutting chunks custom based on position for last 18 regions

# input borders of the regions to be validated

borders = matrix(c(53113091,29288170,73291903,176525230,9795,98395342,80260648,130226803,32603788,81170578,28642588, 81955643,22010347,
                   57376926,33177081,76700977,180783572,2715620,101641945,85190202,134079950,37354675,85001645,33646071,85727849,25354138),  ncol = 2,nrow = 13)
borders = cbind(borders, c(830, 146,534,474,201,1341,1528,758,923,932,1815,136,315))
borders = cbind(borders, c(8, 14,5,4,20,13,15,7,9,9,18,1,3))
colnames(borders) = c("start", "end", "chrchunk", "chr")


for (i in 2:13){

  chr = as.numeric(borders[i, 4])
  chrregion = as.numeric(borders[i, 3])
  
#read position file
  position.file = read.table(paste0("/storage/nipm/kerimbae/swd/swdrecode/chr",chr,".recode.pos.txt"))
  position.file = data.frame(position.file)

#subset position file to a required range
  chunk.pos <- subset(position.file, (position.file$V2 < borders[i,2]& position.file$V2 >borders[i,1])) 
  dim(chunk.pos)[1]
  start = as.numeric(row.names(chunk.pos[1,]))
  end = as.numeric(row.names(chunk.pos[dim(chunk.pos)[1],]))
  length = end-start+1
  print(length)

#load recoded data
 path = paste("/storage/nipm/kerimbae/swd/swdrecode/chr", chr, ".txt",sep ="")
 conn <- file(path, open="r")
 lines <- readLines(conn)

 name <- paste("c", chrregion, sep ="")
 assign(name, lines[start:end])

 write.table(get(name), file = paste("/storage/nipm/kerimbae/validation/chr",chrregion,".txt",sep =""), sep =",",quote = F, row.names = F, col.names = F)
 close(conn)

}