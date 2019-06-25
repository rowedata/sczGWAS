qsub -I -l ncpus=1,mem=2gb,cput=5:0:0 -l walltime=5:0:0 /bin/bash
module load intel intelmpi R
R

#PCA and binary preprocessing

eigenvec <- read.table("/storage/nipm/kerimbae/mgs_scz/qcvcf.eigenvec") 
#5334 by 12 (don't know what first 2 columns mean)

# extract sex data
fam <- read.table("/storage/nipm/kerimbae/mgs_scz/scz.fam")
#5334 by 6 Family ID Sample ID Paternal ID Maternal ID Sex (1=male; 2=female; other=unknown) Affection (0=unknown; 1=unaffected; 2=affected)

#recode sex to 0/1 (0-male)
 
sex <- fam$V5
sex<-replace(sex,sex=='1',0)
sex<-replace(sex,sex=='2',1)

#compose a design matrix

design <-cbind(eigenvec[,3:12],sex)
write.table(design, "/storage/nipm/kerimbae/mgs_scz/design")

#isolate the response vector
  
#recode response to 0/1 (1-affected)  
response <- fam$V6

response <- replace(response, response=='1',0)
response <- replace(response, response=='2',1)

write.table(response,"/storage/nipm/kerimbae/mgs_scz/response")

#############################################################################
# part 2 fit logistic regression and analyze errors. 

diagnosis <- read.table("/storage/nipm/kerimbae/mgs_scz/response",stringsAsFactors=FALSE)
diagnosis <- matrix(unlist(diagnosis), ncol = 1, byrow = TRUE)

covariates <- read.table("/storage/nipm/kerimbae/mgs_scz/design",stringsAsFactors=FALSE)
covariates <- matrix(unlist(covariates), ncol = 11, byrow =FALSE)
#write.table(covariates,"/storage/nipm/kerimbae/mgs_scz/covariates")

#fitting logit

logitmod <- glm(diagnosis ~ covariates, family=binomial(link=logit))
summary(logitmod)
pchisq(deviance(logitmod),df.residual(logitmod),lower=FALSE) #good fit if >0.05
confint(logitmod)

#fitting probit
library(MASS)
probitmod <- glm(diagnosis ~ covariates, family=binomial(link=probit))

summary(probitmod)
pchisq(deviance(probitmod),df.residual(probitmod),lower=FALSE) #good fit if >0.05
confint(probitmod)

#prediction
# We show how to predict the response at dose of 2.5:
x0 <- c(1,covariates[1,])
eta0 <- sum(covariates*coef(probitmod))
ilogit(eta0)

#########################Working vector Y tilde
#alpha coefficients
alphas <- c(0.2397054,6.0258206,0.9150710,-3.2395109,2.3185416,1.8713296,1.1450852,6.7025039,9.0176767,1.9103886,-0.9178600,-0.5622860)
newdata<- rep(1,5334)
newcov <- cbind(newdata,covariates)
  
#
q <- newcov*alphas #Xi\alpha 
xiAlpha <- rowSums(q)
pi0<- exp(xiAlpha)/1+exp(xiAlpha)

newphenoscz <- (diagnosis-pi0)/pi0*(1-pi0)



write.table(newphenoscz,row.names = FALSE,col.names = FALSE, "/storage/nipm/kerimbae/pimass/input/mgsInput/newphenoscz.txt")








