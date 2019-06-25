#Benazir Rowe Spring 2018 UNLV 
#LD tables production from piMASS
#Takes lists of significant loci available 
#Reads the corresponding full mcmc output files and pick top 10% SNPs
#Save names and locations to a file under a region name

#bash command to start the R session
  
qsub -I -l ncpus=1,mem=5gb,cput=5:0:0 -l walltime=5:0:0 /bin/bash
module load intel intelmpi R
R

#specify the chromosome number to analyze

chr = 9


####################### same code for all chromosomes

table = read.table(paste("routput/list",chr,".txt",sep=""), head = T)
length(table$chunk)
list = unique(table$chunk)  #leave unique elements only
length(list)
list =sort(list)
list



for (i in 1 : length(list))
{
  chunk = list[i]
  maxregion = read.table(paste("progs/pimass/pimass/output/chr",chr,"/goldstand/c",chunk,".mcmc.txt",sep=""),head=T)
  #find max only
  max(maxregion$postc)
  which.max(maxregion$postc)
  maxregion$rs[which.max(maxregion$postc)]
  #find top 10%
  cutoff = quantile(maxregion$postc,probs = 0.99)
  topten = vector()
  
    for (j in 1:1000){
      
        if (maxregion$postc[j] > cutoff) {
          topten = c(topten, j) #indices of values in top 10
        }
    }
  
  ldmaterial = cbind(chr,chunk, as.character(maxregion$rs[topten]),maxregion$pos[topten]) #gives a list of snps in top 10%
  class(ldmaterial) #factor
  
  #now print out and save result, then check Ld with known markers
  
  write.csv(ldmaterial,file = paste("routput/LD",chr,"/LDchr",chr,"chunk",chunk,".csv",sep = ""), row.names = F,quote = F)
  write.table(ldmaterial,file = paste("routput/LDchr",chr,"combo.csv",sep = ""), row.names = FALSE,col.names = FALSE,  quote = F,append=TRUE,sep=",")

}

