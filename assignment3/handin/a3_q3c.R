#Testing the EM algorithm for GMM
source("~/Documents/mcgill/math680/assignment3/mixtureEM.R")

ab = read.csv("~/Documents/mcgill/math680/assignment3/abalone.csv")
colnames(ab) = c("sex","length","diameter","height","wholeweight","shuckedweight",
                 "viscweight","shellweight","rings")

#Taking a random subset to speed things up
set.seed(123)
sampsize = 500
dss = ab[sample(1:dim(ab)[1],sampsize),]

#First algorithm
em_dss = mixtureEM(cbind(dss$length,dss$diameter), C=3, tol=1e-7)

#plot(dss$length, dss$diameter, col=dss$sex)
#labels1=apply(em_dss$P.mat, 1,  which.max)
#points(dss$length, dss$diameter, pch=c("1", "2", "3")[labels1])

#numlab = 1*(dss$sex=="I") + 2*(dss$sex=="F") + 3*(dss$sex=="M")
#sum(numlab==labels1)/length(numlab)

#Second algorithm
em_dss_csig = mixtureEM_csig(cbind(dss$length,dss$diameter), C=3, tol=1e-7)
#plot(dss$length, dss$diameter, col=dss$sex)
#labels2=apply(em_dss_csig$P.mat, 1,  which.max)
#points(dss$length, dss$diameter, pch=c("1", "2", "3")[labels2])

#numlab2 = 2*(dss$sex=="I") + 1*(dss$sex=="F") + 3*(dss$sex=="M")
#sum(numlab2==labels2)/length(numlab2)


