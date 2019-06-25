qsub -I -l ncpus=1,mem=10gb,cput=5:0:0 -l walltime=5:0:0 /bin/bash
module load intel intelmpi R
R


chr15 <-read.table("datafiles/chr15",header = T)
relevant15 <- subset(chr15, (BP>83107801 & BP<86887657))
sorted15 <- relevant15[order(relevant15$P),]

chr19 <- read.table("datafiles/chr19",header=T)
relevant19 <- subset(chr19,(BP>15724023 &BP<22638628))
sorted19 <-relevant19[order(relevant19$P),]


chr9 <- read.table("datafiles/chr9",header=T)
relevant9 <- subset(chr9,(BP>34905605 &BP<70379322))
sorted9 <-relevant9[order(relevant9$P),]


PGC <- read.table("datafiles/SCZvsCONT.sumstats",header = T)

chr8 <- subset(PGC, CHR==8)
relevant8<- subset(chr8,(BP>53113091 &BP<57376926))
sorted8 <-relevant8[order(relevant8$P),]


chr14 <- subset(PGC, CHR==14)

#146
relevant14<- subset(chr14,(BP>29288170 &BP<33177081))
sorted14 <-relevant14[order(relevant14$P),]

#1433
relevant14<- subset(chr14,(BP> 86399092&BP<90573122))
sorted14 <-relevant14[order(relevant14$P),]


chr5 <- subset(PGC, CHR==5)
relevant5<- subset(chr5,(BP>73291903 &BP<76700977))
sorted5 <-relevant5[order(relevant5$P),]
head(sorted5)

chr15 <- subset(PGC, CHR==15)
write.table(chr15, "datafiles/chr15")
chr19 <- subset(PGC, CHR==19)
write.table(chr15, "datafiles/chr19")

#analysis october 15
chr9 <- subset(PGC, CHR==9)
c1<- subset(chr9,(BP>34905605 &BP<70379322))#24
c2<- subset(chr9,(BP>32603788 &BP<37354675))#23
c3<- subset(chr9,(BP>81170578 &BP<85001645))#32

chr15<-subset(PGC, CHR==15)

c1<- subset(chr15,(BP>83907801 &BP<86887657))#29
c2<- subset(chr15,(BP>80260648 &BP<85190202)) #28

chr19 <-subset(PGC, CHR==19)
c1<- subset(chr19,(BP>15724023 &BP<22638628))#5

chr8 <-subset(PGC, CHR==8)
c1<- subset(chr8,(BP>53113091 &BP<57376926))#30

chr14<-subset(PGC, CHR==14)

c1<- subset(chr14,(BP>86399092 &BP<90573122))#33
c2<- subset(chr14,(BP>29288170 &BP<33177081)) #6

chr5 <-subset(PGC, CHR==5)
c1<- subset(chr5,(BP>73291903 &BP<76700977))#34


chr4 <-subset(PGC, CHR==4)
c1<- subset(chr4,(BP>176525230 &BP<180783572))#74

chr20 <-subset(PGC, CHR==20)
c1<- subset(chr20,(BP>9795 &BP<2715620))#1


chr13 <-subset(PGC, CHR==13)
c1<- subset(chr13,(BP>98395342 &BP<101641945))#41

chr7 <-subset(PGC, CHR==7)
c1<- subset(chr7,(BP>130226803 &BP<134079950))#58

chr18 <-subset(PGC, CHR==18)
c1<- subset(chr18,(BP>28642588 &BP<33646071))#15

chr1 <-subset(PGC, CHR==1)
c1<- subset(chr1,(BP>81955643 &BP<85727849))#36

chr3 <-subset(PGC, CHR==3)
c1<- subset(chr3,(BP>22010347 &BP<25354138))#15

chr16 <-subset(PGC, CHR==16)
c1<- subset(chr16,(BP> &BP<5512127))#16

write.table(chr15, "datafiles/chr9")

