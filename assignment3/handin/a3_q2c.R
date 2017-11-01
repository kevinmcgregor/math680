#Testing bridgeLogistic function
library(boot)

set.seed(500)
n=50
p=75

X=cbind(1, matrix(rnorm(n*(p-1)), nrow=n, ncol=(p-1)))
beta.star=rnorm(p, mean=0, sd=1/sqrt(p))
n.list=rep(30, n)

#Success prob
pi.list=inv.logit(as.numeric(X%*%beta.star))
#Response
y = rbinom(n=n, size=n.list, prob=pi.list)/n.list

bl = bridgeLogistic(X,y,n.list)

