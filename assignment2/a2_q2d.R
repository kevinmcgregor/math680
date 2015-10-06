source("~/Documents/mcgill/math680/assignment2/centeredRidge_func.R")

#Testing out the function (borrowed from Notes 2)
set.seed(680)
n=50; p=10; lambda=1; sigma=1
X=cbind(1, matrix(rnorm(n*p), nrow=n, ncol=p))
Xtilde = scale(X[,-1], scale=FALSE)
beta.star=c(1,runif(p,-5,5))
y=X%*%beta.star + sigma*rnorm(n)
centeredRidge(X,y,lambda)

#Trying replications and comparing to theoretical variance of estimator
nrep = 1000
betarep = matrix(0,p,nrep)
ycur = NULL
for (i in 1:nrep) {
  ycur=X%*%beta.star + sigma*rnorm(n)
  betarep[,i] = centeredRidge(X,ycur,lambda)$b
}

betaexpect = apply(betarep,1,mean)
theor_betaexpect = qr.solve(crossprod(Xtilde, Xtilde)+lambda*diag(p))%*%crossprod(Xtilde, Xtilde)%*%beta.star[-1]

betavar = (n-1)/n * cov(t(betarep))
theor_betavar = sigma*qr.solve(crossprod(Xtilde, Xtilde)+lambda*diag(p))%*%t(Xtilde)%*%
                    t(qr.solve(crossprod(Xtilde, Xtilde)+lambda*diag(p))%*%t(Xtilde))
