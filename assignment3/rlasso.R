###################################################################
## Compute the lasso penalized least squares solution
###################################################################
## arguments: 
##  X, the n by p design matrix with all entries in
##     its first column equal to one.
##  y, the response vector with n entries.
##  lam, the non-negative tuning parameter value (or a vector of p values)
##  tol, convergence tolerance
##  maxit, the maximum number of iterations allowed
##  maxit.in, the maximum number of inner loop iterations allowed
##  quiet, should the function stay quiet?
##  b, optional 0th iterate for the regression 
##           coefficient estimates
##  h, optional 0th iterate value for the h vector 
##           (required if b.start is unspecified)
##  fval, optional 0th iterate objective function value 
##              (required if b.start is unspecified)
##  XtX, optional, useful when computing the solution at multple lambdas 
##  Xty, optional, useful when computing the solution at multple lambdas 
##  xbar, optional, useful when computing the solution at multple lambdas 
##  ybar, optional, useful when computing the solution at multple lambdas 
##
## this function returns a list with nine elements:
##  b1, is the intercept estimate.
##  b,  is the vector of remaining regression coefficient
##      estimates.  It has p-1 entries.
##  h, this is the h vector at the final iterate
##  fval, this is the objective function value at the final iterate
##  total.iterations, number of iterations made
##  XtX 
##  Xty 
##  xbar 
##  ybar 
###################################################################

minlassoLSwithC=function(X, y, lam, tol=1e-7, maxit=1e4, quiet=FALSE, 
                         b=NULL, h=NULL, fval=NULL,
                         XtX=NULL, Xty=NULL, xbar=NULL, ybar=NULL)
{
 
  p=ncol(X)
  informed=(!is.null(XtX)) & (!is.null(Xty)) & (!is.null(xbar)) & (!is.null(ybar))
  if(!informed)
  {
    ## compute the column-centered predictor matrix Xc
    xbar=apply(X[,-1, drop=FALSE], 2, mean)
    Xc=scale(X[,-1, drop=FALSE], scale=FALSE, center=xbar)
  
    ## compute the centered response vector yc
    ybar=mean(y)
    yc=y-ybar
  
    ## compute matrices for later use
    XtX=crossprod(Xc)
    Xty=crossprod(Xc,yc)
  }
  
  ## use zeros starting value if unspecified
  warm.start=(!is.null(b)) & (!is.null(h)) & (!is.null(fval))
  if(!warm.start)
  {
    b=rep(0, p-1)
    h=rep(0, p-1)
    fval=0.5*sum(yc^2)
  } 
  
  tol.adjusted=fval*tol  
  pminus1=p-1
  
  if(!is.loaded("rlasso")) dyn.load("rlasso.so")
  total=0
  dotCoutput=.C("rlasso", b=as.double(b), Xty=as.double(Xty), XtX=as.double(XtX), 
                          h=as.double(h), pin=as.integer(pminus1), lam=as.double(lam), 
                          tol=as.double(tol.adjusted), maxit=as.integer(maxit), 
                          objective=as.double(fval), total=as.integer(total))
  if(!quiet)
    cat("It took", dotCoutput$total,  "iterations\n")
  
  ## compute intercept estimate
  b1=ybar - sum(xbar * dotCoutput$b)
  
  return(list(b1=b1, b=dotCoutput$b, h=dotCoutput$h, 
              fval=dotCoutput$objective, 
              total.iterations=dotCoutput$total, 
              XtX=XtX, Xty=Xty, ybar=ybar, xbar=xbar))
} 

 
 


soft=function(u, lam)
{
  tmp=abs(u) - lam
  tmp = tmp*(tmp > 0)
  val=sign(u)*tmp
  return(val)
}
###################################################################
## The R version of the function above
###################################################################
minlassoLS=function(X, y, lam, tol=1e-7, maxit=1e4, quiet=FALSE, 
                         b=NULL, h=NULL, fval=NULL,
                         XtX=NULL, Xty=NULL, xbar=NULL, ybar=NULL)
{
 
  p=ncol(X)
  informed=(!is.null(XtX)) & (!is.null(Xty)) & (!is.null(xbar)) & (!is.null(ybar))
  if(!informed)
  {
    ## compute the column-centered predictor matrix Xc
    xbar=apply(X[,-1, drop=FALSE], 2, mean)
    Xc=scale(X[,-1, drop=FALSE], scale=FALSE, center=xbar)
  
    ## compute the centered response vector yc
    ybar=mean(y)
    yc=y-ybar
  
    ## compute matrices for later use
    XtX=crossprod(Xc)
    Xty=crossprod(Xc,yc)
  }
  
  ## use zeros starting value if unspecified
  warm.start=(!is.null(b)) & (!is.null(h)) & (!is.null(fval))
  if(!warm.start)
  {
    b=rep(0, p-1)
    h=rep(0, p-1)
    fval=0.5*sum(yc^2)
  } 
  
  iterating=TRUE
  k=0
  tol.adjusted=fval*tol  
  if(!quiet) cat("k=", k, "fval=", fval, "\n") 
    
  while(iterating)
  {
    k=k+1
    
    ## total.diff will be the total absolute difference
    ## from advancing all (p-1) entries in our iterate
   
    totalchange=0
   
    for( j in 1:(p-1) )
    {
      ## update the jth entry
      bj.new=soft(u=(Xty[j] - h[j]), lam=lam)/XtX[j,j]   
      if(bj.new != b[j] )
      {
        i=1:(p-1)
        i=i[-j]
        
        ## compute useful differences
        db=b[j] - bj.new
        dh=-XtX[i,j]*db
        
        ## compute objective function value change 
        dfval=db*(Xty[j] -0.5*h[j]) + 0.5*XtX[j,j]*(bj.new^2 - b[j]^2)
        dfval=dfval + 0.5*sum(b[i]*dh)+lam*( abs(bj.new)-abs(b[j]) )
        
        ## update current objective function value
        fval=fval+dfval
        
        ## update total change
        totalchange=totalchange-dfval
        
        ## update b[j]
        b[j]=bj.new
        
        ## update the vector h[i]
        h[i]= h[i] + dh
      }
    }
    
    if( ( totalchange < tol.adjusted) | (k >= maxit)) 
      iterating=FALSE

    if(!quiet) cat("k=", k, "fval=", fval, "\n") 
  }
  
  ## compute intercept estimate
  b1=ybar - sum(xbar * b)
  
  ## make b a vector instead of a one column matrix
  b=as.numeric(b)
  
  return(list(b1=b1, b=b, h=h, fval=fval, total.iterations=k,
              XtX=XtX, Xty=Xty, ybar=ybar, xbar=xbar))
} 

## use the byte code compiler
library(compiler)
minlassoLSbcc=cmpfun(minlassoLS)


