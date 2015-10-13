#Testing out the normal ridge regression function

#Testing out the function (borrowed from Notes 2)
set.seed(680)
n=100; p=10; lambda=1; sigma=1
X=cbind(1, matrix(rnorm(n*p), nrow=n, ncol=p))
Xtilde = scale(X[,-1], scale=FALSE)
beta.star=c(1,runif(p,-5,5))
y=X%*%beta.star + sigma*rnorm(n)
ytilde = scale(y, scale=FALSE)

#Testing the normal ridge regression function
tmp = normRidge(X,y,lambda)

#Checking gradient
beta.nr = tmp$b
sigma.nr = tmp$sigma2
grad.beta = (crossprod(Xtilde)%*%beta.nr - crossprod(Xtilde, ytilde))/sigma.nr + lambda*beta.nr
grad.sigma = n/sigma.nr - (crossprod(ytilde)-2*crossprod(Xtilde%*%beta.nr,ytilde)+crossprod(Xtilde%*%beta.nr,Xtilde)%*%beta.nr)/(sigma.nr)^2

