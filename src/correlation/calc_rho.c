#include "grasp.h"
static char *rcsid="Development code";

#define RC(row,col,size) ((row)+(size)*(col))   /* clapack packing conventions  */

typedef float real;                         /* from f2c.h */
typedef struct { real r, i; } complex;      /* from f2c.h */

int calc_rho(int offset,int correlation_width,float threshold,float *rp_signal,float *ip_signal,
	     int nenv_chan,float **rp_env,float **ip_env,float *rho2_pairwise,complex *A,complex *B,float *modx2sum)
{ 
  int k,n,row,col;

  n=nenv_chan; /* Size of matrix */

  /*************************************************************/
  /*                                                           */
  /*   Construct averages over band length correlation_width   */
  /*                                                           */
  /*************************************************************/ 
 
  /* Initialise all sums to zero ... */

  for (col=0;col<n;col++) {
    for (row=0;row<=col;row++) { /* as we only need give the upper triangular part of the matrix */
 	 A[RC(row,col,n)].r=A[RC(row,col,n)].i=0.0;
    }
  }
  for (row=0;row<n;row++) {
      B[row].r=B[row].i=0.0;
    }
  *modx2sum=0.0;

  /* ... and now form sums over bandwidth correlation_width */
  for (k=offset;k<offset+correlation_width;k++) {
    
    for (col=0;col<n;col++) {
	  for (row=0;row<=col;row++) {  /* as we only need give the upper triangular part of the matrix */
	    A[RC(row,col,n)].r+=rp_env[row][k]*rp_env[col][k]+ip_env[row][k]*ip_env[col][k]; /* Re(Y_p^* Y_q) */
	    A[RC(row,col,n)].i+=rp_env[row][k]*ip_env[col][k]-ip_env[row][k]*rp_env[col][k]; /* Im(Y_p^* Y_q) */
	  }
    }
    
    for (row=0;row<n;row++) {
      B[row].r+=rp_signal[k]*rp_env[row][k]+ip_signal[k]*ip_env[row][k];  /* Re(Y_row^* X) */
      B[row].i+=ip_signal[k]*rp_env[row][k]-rp_signal[k]*ip_env[row][k];  /* Im(Y_row^* X) */
    }
    *modx2sum+=rp_signal[k]*rp_signal[k]+ip_signal[k]*ip_signal[k];
  }

  /*  do we threshold */
  for (col=0;col<n;col++) {
    for (row=0;row<=col;row++) {  /* as we only need give the upper triangular part of the matrix */
      
      if ( (B[col].r*B[col].r+B[col].i*B[col].i) < threshold*A[RC(col,col,n)].r*(*modx2sum) ) 
	B[col].r=B[col].i=0.0; 
   
      if ( (A[RC(row,col,n)].r*A[RC(row,col,n)].r+A[RC(row,col,n)].i*A[RC(row,col,n)].i) <
	   threshold*A[RC(row,row,n)].r*A[RC(col,col,n)].r ) 
	A[RC(row,col,n)].r=A[RC(row,col,n)].i=0.0;
    }
        
  }

  
  /* Reduction in variance using the best r for just the single channel (after thresholding) */ 
  for (row=0;row<nenv_chan;row++) {
    rho2_pairwise[row]=(B[row].r*B[row].r+B[row].i*B[row].i)/(A[RC(row,row,n)].r*(*modx2sum));
  }                           

  return(0);
}
