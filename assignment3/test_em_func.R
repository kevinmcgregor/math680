source("~/Documents/mcgill/math680/assignment3/mixtureEM.r")
set.seed(5)
n=1000
p=2
C=3

## create the three covariance matrices
## Sigma 1 is diagonal with equal diagonal elements
Sigma1=0.05*diag(p)

## Sigma 2 is not diagonal:
Evectors=cbind(c(1,1), c(1,-1))/sqrt(2)
Sigma2=Evectors%*% diag(c(0.001, 0.1))%*%t(Evectors)

## Sigma 3 is diagonal with unequal diagonal elements
Sigma3=diag(c(0.1, 0.001))


## compute their matrix square roots for data generation
Sigma1.sqrt=sqrt(0.05)*diag(p)
Sigma2.sqrt=Evectors%*% diag(sqrt(c(0.001, 0.1)))%*%t(Evectors) 
Sigma3.sqrt=diag(sqrt(c(0.1, 0.001)))


## Create the three mu's
mu1=c(0,0)
mu2=c(1,0)
mu3=c(0,1)

## Create the three pi's
pi1=1/3
pi2=1/3
pi3=1/3


## Generate the data matrix:
X=matrix(NA, nrow=n, ncol=p)
y=rep(0,n)
for(i in 1:n)
{
  ## perform a multinomial trial to
  ## generate the response cateogory
  mtrial=rmultinom(1, size=1, prob=c(pi1,pi2,pi3))
  if(mtrial[1])  ## resulted in category 1
  {
    X[i,] = mu1 + Sigma1.sqrt%*%rnorm(p)
    y[i]=1
  } else if(mtrial[2])  ## resulted in category 2
  {
    X[i,] = mu2 + Sigma2.sqrt%*%rnorm(p)
    y[i]=2
  } else  ## resulted in category 3
  {  
    X[i,] = mu3 + Sigma3.sqrt%*%rnorm(p)
    y[i]=3
  }
}

## plot the points without class labels
plot(X)

## plot points with class labels
plot(X, col=y)

## fit the Gaussian mixture model
outfast=mixtureEM(X=X, C=C, tol=1e-7)

## get the assigned cluster/class labels
labels=apply(outfast$P.mat, 1,  which.max)

## add these labels to the plot
points(X, pch=c("1", "2", "3")[labels])

#Trying out the version where sigma is constant over the different categories...
X1=matrix(NA, nrow=n, ncol=p)
y1=rep(0,n)
for(i in 1:n)
{
  ## perform a multinomial trial to
  ## generate the response cateogory
  ## SIGMA EQUAL BETWEEN GROUPS
  mtrial1=rmultinom(1, size=1, prob=c(pi1,pi2,pi3))
  if(mtrial1[1])  ## resulted in category 1
  {
    X1[i,] = mu1 + Sigma2.sqrt%*%rnorm(p)
    y1[i]=1
  } else if(mtrial1[2])  ## resulted in category 2
  {
    X1[i,] = mu2 + Sigma2.sqrt%*%rnorm(p)
    y1[i]=2
  } else  ## resulted in category 3
  {  
    X1[i,] = mu3 + Sigma2.sqrt%*%rnorm(p)
    y1[i]=3
  }
}

plot(X1, col=y1)

outfast_csig=mixtureEM_csig(X=X1, C=C, tol=1e-7)

labels1=apply(outfast_csig$P.mat, 1,  which.max)
points(X1, pch=c("1", "2", "3")[labels1])


