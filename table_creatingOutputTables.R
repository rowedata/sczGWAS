# Benazir Rowe 
# Summer 2018 UNLV 
# Code to produce single true phenotype run summary on Rpubs


# bash command to start the R session on the remote cluster for high volume computations

qsub -I -l ncpus=1,mem=5gb,nmics=1,cput=5:0:0 -l walltime=5:0:0 /bin/bash
module load intel intelmpi R
R


######################################################################################
#change input ONLY here

p = 1000 # number of simulations available
chr = 8  # chromosome number
length = 72  # chromosome length

#########################################
# gold standard

gs <- matrix(,length,5)

  for (i in 1:length){
    
    infile <- paste("progs/pimass/pimass/output/chr",chr,"/goldstand/c",i,".mcmc.txt", sep = "")
    mcmc <- read.table(infile, head = T, stringsAsFactors = FALSE)
    
    gs[i,1] <- sum(mcmc$postc)
    gs[i,2] <- mean(mcmc$postc)
    gs[i,3] <- max(mcmc$postc)
    gs[i,4] <- min(mcmc$pos)
    gs[i,5] <- max(mcmc$pos)
  }

colnames(gs) = c("sum","mean","max", "start","end")
write.table(gs, file =paste("routput/gschr",chr,".txt",sep=""), row.names = F, col.names = T)
gs


##########################################################################################
#analyzing 1000 simulations

for (k in 1:length){
  
  sim <- matrix(, p, 3) #value
  colnames(sim) <- c("sum", "mean", "max")
  
  for (i in 1:p){
    infile <- paste("progs/pimass/pimass/output/chr",chr,"/p",i,"/c", k,".mcmc.txt", sep = "")

    mcmc <- read.table(infile, head = T, stringsAsFactors = FALSE)
    sim[i, 1] <- sum(mcmc$postc, na.rm = TRUE)
    sim[i, 2] <- mean(mcmc$postc,na.rm = TRUE)
    sim[i, 3] <- max(mcmc$postc,na.rm = TRUE)
  }
  
  
  j <- noquote(paste("sim", k ,sep ="",collapse = NULL))  #name
  assign(j, sim)                                          #assign value to name
}

#####################################################
pval <- matrix(, length, 5)

for (k in 1:length){
  
  w <- noquote(paste("sim", k ,sep ="",collapse = NULL))
  sim <-   get(w)
  
  for (j in 1:3){

    sum(sim[, j] > gs[k, j])
    pval[k, j] <- sum(sim[, j] > gs[k, j]) / p
  }
  pval[k, 4] = round(gs[k,4]/1000000, 1) #start point in megabases
  pval[k, 5] = round(gs[k,5]/1000000, 1) #end point in megabases
  
}

colnames(pval) <- c("sum","mean","max", "start","end")
write.table(pval,file = paste("routput/pval", chr, ".txt",sep=""), row.names = F, col.names = T)


#create table with regions that exceed alpha cutoff

alpha = 0.05

x = vector()
xi = vector() 
xj = vector() 

for (j in 1:3){
  for (i in 1:length){
    if (pval[i,j] < alpha) {
      x = c(x,colnames(pval )[j])
      xi = c(xi,i)
      xj = c(xj,pval[i,j])
    }
  }
}

q <- cbind(x, xi, xj)
colnames(q) <- c("metric","chunk", "pvalue")
write.table(q, file = paste("routput/list", chr, ".txt",sep=""), row.names = F, col.names = T)



