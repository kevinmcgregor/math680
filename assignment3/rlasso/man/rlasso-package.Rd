\name{rlasso-package}
\alias{rlasso-package}
\alias{rlasso}
\docType{package}
\title{
Functions for doing LASSO
}
\description{
Functions for doing LASSO
}
\details{

\code{minlassoLS} performs LASSO-penalized regression.
\code{minlassoLSwithC} performs LASSO-penalized regression using an external C function.
}
\author{
Yi Yang, Kevin McGregor
Maintainer: Kevin McGregor
}
\references{
}
\keyword{ LASSO }
\seealso{
}
\examples{
  n=50
  X = cbind(1,rnorm(n,1,2),rnorm(5,1))
  beta.star = rep(1,3)
  y = X\%*\%beta.star + rnorm(n,0,1)
  
  fit = minlassoLS(X, y, lam=1)
  fit2 = minlassoLSwithC(X, y, lam=1)

}
