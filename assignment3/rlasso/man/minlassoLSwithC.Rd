\name{minlassoLSwithC}
\alias{minlassoLSwithC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Compute the lasso penalized least squares solution using C
}
\description{
Compute the lasso penalized least squares solution using C
}
\usage{
minlassoLSwithC(X, y, lam, tol = 1e-07, maxit = 10000, quiet = FALSE, b = NULL, h = NULL, fval = NULL, XtX = NULL, Xty = NULL, xbar = NULL, ybar = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
  The n by p design matrix with all entries in its first column equal to one.
}
  \item{y}{
  The response vector with n entries.
}
  \item{lam}{
  The non-negative tuning parameter value (or a vector of p values).
}
  \item{tol}{
  Convergence tolerance.
}
  \item{maxit}{
  The maximum number of iterations allowed.
}
  \item{quiet}{
  Should the function stay quiet?
}
  \item{b}{
  optional 0th iterate for the regression coefficient estimates.
}
  \item{h}{
  optional 0th iterate value for the h vector (required if b.start is unspecified)
}
  \item{fval}{
  optional 0th iterate objective function value (required if b.start is unspecified)
}
  \item{XtX}{
  optional, useful when computing the solution at multple lambdas
}
  \item{Xty}{
  optional, useful when computing the solution at multple lambdas 
}
  \item{xbar}{
  optional, useful when computing the solution at multple lambdas
}
  \item{ybar}{
  optional, useful when computing the solution at multple lambdas
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{b1}{ the intercept estimate}
  \item{b}{is the vector of remaining regression coefficient estimates.  It has p-1 entries.}
  \item{h}{this is the h vector at the final iterate}
  \item{fval}{this is the objective function value at the final iterate}
  \item{total.iterations}{number of iterations made}
  \item{XtX}{} 
  \item{Xty}{} 
  \item{xbar}{} 
  \item{ybar}{}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Yi Yang, Kevin McGregor
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
  n=50
  X = cbind(1,rnorm(n,1,2),rnorm(5,1))
  beta.star = rep(1,3)
  y = X\%*\%beta.star + rnorm(n,0,1)
  fit2 = minlassoLSwithC(X, y, lam=1)
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ C }
\keyword{ LASSO }% __ONLY ONE__ keyword per line
