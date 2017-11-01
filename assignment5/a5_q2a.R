#Testing the fused ridge cross validation function (cvfr)

source("~/Documents/mcgill/math680/assignment5/fusedRidge_crossval.R")

set.seed(680)
n=50; p=10; sigma=1
X=cbind(1, matrix(rnorm(n*p), nrow=n, ncol=p))
Xtilde = scale(X[,-1], scale=FALSE)
beta.star=c(1,runif(p,-5,5))
y=X%*%beta.star + sigma*rnorm(n)

lam = 10^(seq(-8,8, by=0.5))

f_fit = cvfr(X, y, lam, lam, 2)




