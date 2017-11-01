#Testing out the bridge regression function
source("~/Documents/mcgill/math680/assignment3/bridgeReg.R")

#Testing out the function (data generation borrowed from notes)
set.seed(123)
n=200
p=50

omega.star=2
## generate the true beta (so some of its entries are zero)
#beta.star=rbinom(p, size=1, prob=0.3)*1
beta.star = rnorm(p,0,2)
th=0.3
sigma=diag(p-1)
for(i in 1:(p-1)) for(j in 1:(p-1)) 
  sigma[i,j]=th^abs(i-j)
eo=eigen(sigma, symmetric=TRUE)
sigma.sqrt=eo$vec%*%diag(eo$val^0.5)%*%t(eo$vec)
predictor.matrix=matrix(rnorm(n*(p-1)), nrow=n, ncol=(p-1))%*%sigma.sqrt
X=cbind(1, predictor.matrix) 
Xc = scale(X[,-1], scale=FALSE)
y=X%*%beta.star + rnorm(n, mean=0, sd=sqrt(1/omega.star))
yc = scale(y, scale=FALSE)

lambda=1
alpha=1.5
br = bridgeReg(y,X,lam=1, alpha=1.5)
beta.br = br$b

grad.br = -crossprod(Xc,y) + crossprod(Xc)%*%beta.br + lambda*sign(beta.br)*abs(beta.br)^(alpha-1)

