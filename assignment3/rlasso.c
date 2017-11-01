#include <math.h>
 
/*
To compile 

MAKEFLAGS="CFLAGS=-fPIC" 
gcc -c rlasso.c -fPIC
R CMD SHLIB rlasso.o
 
*/

void rlasso(double * b, double * Xty, double *XtX, double *h, int * pin, 
            double * lam, double * tol, int * maxit, double * objective, 
	    int * total)
{
  double fval=*objective;
  int p=*pin;
  int k, j, jp, jppj, i;
  double bnew, bold, dh, db, dfval, soft_arg, 
         abs_soft_arg, tmp, sum_b_dh, total_diff;
  
  k=0;
  total_diff=*tol+1.0;
  
  while( (total_diff > *tol) && (k < *maxit) )
  {
    k++;
    total_diff = 0.0;
    
    for(j=0; j <p; j++)
    {
      //use these to avoid redundant flops
      jp=j*p;
      jppj=jp+j;
      
      /* compute the update for b[j]
         with the soft thresholding operation */  
      soft_arg=Xty[j]-h[j];
      abs_soft_arg=fabs(soft_arg);
      tmp=abs_soft_arg - *lam;
      bnew=0.0;
      if(tmp  > 0.0 )
      {
        if(soft_arg > 0.0)
          bnew = tmp;
        else if( soft_arg < 0.0 )
          bnew = -tmp;
        else
          bnew=0.0;
      }            
      bnew=bnew/XtX[jppj];
      
      
      if(bnew != b[j])
      {
	// compute useful quantities
	bold=b[j];
	db=bold-bnew;
	
	// update h and compute a sum used to update fval
	sum_b_dh=0.0;
	for(i=0; i < p; i++)
	{
	  if(i != j)
	  {  
	    dh=-XtX[jp+i]*db;
	    // update the ith element of h
	    h[i]+=dh;
	    // add to the sum
	    sum_b_dh+=b[i]*dh;
	  }
	} 
	
	// compute objective function value change 
        dfval=db*(Xty[j] -0.5*h[j]) + 0.5*XtX[jppj]*(bnew*bnew - bold*bold);
        dfval+= 0.5*sum_b_dh + (*lam)*( fabs(bnew)-fabs(bold) );
	
	// update the objective function value
        fval+=dfval;
	
	// update how much it has changed in this iteration
        total_diff-=dfval;
	
	// update b[j]
	b[j]=bnew;
	
      } // end if 
      
    } // end for loop over j
    
  } // end while loop
  
  // send information back to R:
  total[0]=k;
  objective[0]=fval;
}
 
