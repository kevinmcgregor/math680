#Testing R function for A4 question 1

source("fitBetaNeg.R")

#Testing basic case
set.seed(234908)
n=50
p=4
X = cbind(1,rnorm(n),rnorm(n,4,2),rnorm(n,3,3))
beta.star = c(1,1,1,6)
y = X%*%beta.star + rnorm(n)

fit1 = fitBetaNeg(X,y)

#part (b) Testing normal cases model when beta.star_p is positive

set.seed(1203)
n=50
p=5

omega.star=2
beta.star = runif(p,0,5)
th=0.3
sigma=diag(p-1)
for(i in 1:(p-1)) for(j in 1:(p-1)) 
  sigma[i,j]=th^abs(i-j)
eo=eigen(sigma, symmetric=TRUE)
sigma.sqrt=eo$vec%*%diag(eo$val^0.5)%*%t(eo$vec)
predictor.matrix=matrix(rnorm(n*(p-1)), nrow=n, ncol=(p-1))%*%sigma.sqrt
X=cbind(1, predictor.matrix) 
y=X%*%beta.star + rnorm(n, mean=0, sd=sqrt(1/omega.star))

fit_partb = fitBetaNeg(X,y)
#Checking KKT conditions
R=c(rep(0,p-1),1)
u=2*t(R)%*%qr.solve(crossprod(X))%*%crossprod(X,y)/(crossprod(R,qr.solve(crossprod(X)))%*%R)
-crossprod(X,y) + crossprod(X)%*%fit_partb + u/2*R
crossprod(R, fit_partb)
u*crossprod(R,fit_partb)
u

#part (c) Testing normal cases model when beta.star_p is negative

set.seed(1203)
n=50
p=5

omega.star=2
beta.star = c(runif(p-1,0,4),-2)
th=0.3
sigma=diag(p-1)
for(i in 1:(p-1)) for(j in 1:(p-1)) 
  sigma[i,j]=th^abs(i-j)
eo=eigen(sigma, symmetric=TRUE)
sigma.sqrt=eo$vec%*%diag(eo$val^0.5)%*%t(eo$vec)
predictor.matrix=matrix(rnorm(n*(p-1)), nrow=n, ncol=(p-1))%*%sigma.sqrt
X=cbind(1, predictor.matrix) 
y=X%*%beta.star + rnorm(n, mean=0, sd=sqrt(1/omega.star))

fit_partc = fitBetaNeg(X,y)
#Checking KKT conditions
R=c(rep(0,p-1),1)
u=0 # u is zero since the OLS estimate is already negative
-crossprod(X,y) + crossprod(X)%*%fit_partc + u/2*R
crossprod(R, fit_partc)
u*crossprod(R,fit_partc)
u

