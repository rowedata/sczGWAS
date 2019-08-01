# Benazir Rowe 
# Summer 2018 UNLV 
# Validation region creation based on position of the original MGS dataset
# based on genetic position find the corresponding SNPs in the validation SSCCS dataset

#input the regions start and end positions

borders = matrix(c(34905605,86399092, 83907801,15724023,70379322,90573122,86887657,22638628),  ncol = 2, nrow = 4)
borders = cbind(borders, c(9,14,15,19))
colnames(borders) = c("start", "end", "chr")

i = 4
chr = as.numeric(borders[i,3])

position.file = read.table(paste0("/storage/nipm/kerimbae/swd/swdrecode/chr",chr,".recode.pos.txt"))
position.file = data.frame(position.file)

chunk.pos <- subset(position.file, (position.file$V2 < borders[i, 2]& position.file$V2 > borders[i,1])) 
dim(chunk.pos)[1]

start = as.numeric(row.names(chunk.pos[1, ]))
end = as.numeric(row.names(chunk.pos[dim(chunk.pos)[1],]))
length = end-start+1
print(length)

path = paste("/storage/nipm/kerimbae/swd/swdrecode/chr", chr, ".txt",sep ="") # location to save the file
conn <- file(path,open="r") # open connection to the chosen location
lines <- readLines(conn)    # 

name <- paste("c", chr, sep = "")
assign(name, lines[start:end])

# write.table(get(name), file = paste("/storage/nipm/kerimbae/validation/chr",chr,".txt",sep =""), sep =",",quote = F, row.names = F, col.names = F)
close(conn)


#next 3 regions to validate
#####################################################################

borders = matrix(c(53113091,29288170,73291903,57376926,33177081,76700977),  ncol = 2,nrow = 3)
borders = cbind(borders, c(8,14,5))
borders = cbind(borders,c(30,6,34))
colnames(borders)=c("start","end", "chr","chunk")

i = 3

#retrieve chromosome & chunk info
chr = as.numeric(borders[i,3])
chunk = as.numeric(borders[i,4])

#read position file
position.file = read.table(paste0("/storage/nipm/kerimbae/swd/swdrecode/chr",chr,".recode.pos.txt"))
position.file = data.frame(position.file)
dim(position.file)

#subset position file to a required range
chunk.pos <- subset(position.file, (position.file$V2 < borders[i,2] & position.file$V2 > borders[i,1])) 
dim(chunk.pos)[1]

start = as.numeric(row.names(chunk.pos[1,])) #line number of first cut dataset
end = as.numeric(row.names(chunk.pos[dim(chunk.pos)[1],])) #line number of last cut dataset
length = end - start + 1
print(length)

path = paste("/storage/nipm/kerimbae/swd/swdrecode/chr", chr, ".txt",sep ="")
conn <- file(path,open="r")
lines <- readLines(conn)

name <- paste0("c", chr)
assign(name, lines[start:end])

write.table(get(name), file = paste("/storage/nipm/kerimbae/validation/chr",chr,chunk,".txt",sep =""), sep =",",quote = F, row.names = F, col.names = F)
close(conn)




