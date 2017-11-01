setwd("~/Documents/mcgill/math680/assignment3/")

source("rlasso.r")

set.seed(1)
n=100
p=8

## Define the true parameter values
## the error precision
omega.star=2

## generate the true beta (so some of its entries are zero)
beta.star=rbinom(p, size=1, prob=0.3)*1

## the covariance matrix for the random predictors
## is AR(1), 0.3
th=0.3
sigma=diag(p-1)
for(i in 1:(p-1)) for(j in 1:(p-1)) 
  sigma[i,j]=th^abs(i-j)
  
## compute the square root of the covariance 
## matrix for the random predictors  
eo=eigen(sigma, symmetric=TRUE)
sigma.sqrt=eo$vec%*%diag(eo$val^0.5)%*%t(eo$vec)

## generate a realization of the random predictor matrix
## the rows are independent draws from N_{p-1}(0, sigma)
predictor.matrix=matrix(rnorm(n*(p-1)), nrow=n, ncol=(p-1))%*%sigma.sqrt

## create the design matrix
X=cbind(1, predictor.matrix) 

## generate the response vector
y=X%*%beta.star + rnorm(n, mean=0, sd=sqrt(1/omega.star))

#########################################################
##  Compare three functions to compute the solution:
##  pick a positive value for lam
lam=0.5

## time using the R version
system.time(expr=(fit=minlassoLS(X=X, y=y, lam=lam, tol=1e-10, quiet=TRUE, maxit=1e4)))

## time using the byte-code compiled R version
system.time(expr=(fitbcc=minlassoLSbcc(X=X, y=y, lam=lam, tol=1e-10, quiet=TRUE, maxit=1e4)))

## time using the C version
system.time(expr=(fitc=minlassoLSwithC(X=X, y=y, lam=lam, tol=1e-10, quiet=TRUE, maxit=1e4)))


## check if there is an approximate zero subgradient at the final iterate
xbar=apply(X[,-1, drop=FALSE], 2, mean)
Xc=scale(X[,-1, drop=FALSE], scale=FALSE, center=xbar)
ybar=mean(y)
yc=y-ybar

eta = sign(fitc$b)*(fitc$b != 0) + (1/lam)*(crossprod(Xc, yc)-crossprod(Xc)%*%fitc$b)*(fitc$b == 0)
## the following must be less than or equal to one
max(abs(eta))
subgradatsol=-crossprod(Xc, yc) + crossprod(Xc)%*%fitc$b + lam * eta
## the following should be close to zero
max(abs(subgradatsol))


##########################################################
## compute the solution path
####

lam.vec=10^seq(from=-6, to=log10(max(abs(crossprod(Xc,yc)))), length.out=50)
## use decreasing order:
lam.vec=rev(lam.vec) 

bhat.matrix=matrix(NA, nrow=p, ncol=length(lam.vec))
for(j in 1:length(lam.vec))
{
  if(j==1)
    fitc=NULL
  ## use warm starting values when j > 1  
  fitc=minlassoLSwithC(X=X, y=y, lam=lam.vec[j], b=fitc$b, 
                       h=fitc$h, fval=fitc$fval, 
                       tol=1e-12, quiet=TRUE, maxit=1e4,
                       XtX=fitc$XtX, Xty=fitc$Xty, xbar=fitc$xbar,
                       ybar=fitc$ybar)
  bhat.matrix[,j]=c(fitc$b1, fitc$b)
}

## plot of the estimated coefficients vs 
## the tuning parameter value (exclude the intercept)
matplot(log10(lam.vec), t(bhat.matrix[-1,]), t="o")

###########################################################









