set.seed(680)
n=25;p=10
X=cbind(1, matrix(rnorm(n*(p-1)), nrow=n, ncol=(p-1)))
beta.star=rep(1,p)
y=X%*%beta.star + 1*rnorm(n)
## compute the ols estimate the old way
beta.hat=lm.fit(x=X, y=y)$coeff
## compute the ols estimate the new way
xbar=apply(X[,-1], 2, mean); ybar=mean(y)
Xtilde=scale(X[,-1], center=xbar, scale=FALSE); ytilde=y-ybar
beta.hat.minus1=lm.fit(x=Xtilde, y=ytilde)$coeff
beta.hat.1=ybar - sum(xbar * beta.hat.minus1)
beta.hat.new=c(beta.hat.1, beta.hat.minus1)

hattilde = Xtilde%*%solve(t(Xtilde)%*%Xtilde)%*%t(Xtilde)

hattilde%*%y
fit=lm.fit(x=X, y=y)

(diag(n)-hattilde)%*%rep(1,n)*ybar + hattilde%*%y

rep(1,n)*ybar + hattilde%*%(y - rep(1,n)*ybar)

Xni = X[,-1]
Xnihat = Xni%*%qr.solve(t(Xni)%*%Xni)%*%t(Xni)
Xnihat%*%y


