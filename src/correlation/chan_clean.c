#include "grasp.h"
#include  "f2c.h"
static char *rcsid="Development code";

#define RC(row,col,size) ((row)+(size)*(col))
  
int chan_clean(int offset,int correlation_width,float threshold,float *rho2,float *rp_signal,float *ip_signal,
	       int nenv_chan,float **rp_env,float **ip_env,float *rp_clean,float *ip_clean,complex *A,complex *B,
	       float modx2sum,complex *R,complex *work,integer lwork,integer *ipivot)
{ 
  /*******************************************************************************/
  /*          NB f2c.h defines a data types complex and integer (= long int)     */
  /*                   which is used by clapack routines                         */
  /*******************************************************************************/
  int chesv_(char *, integer *, integer *, complex *,integer *,integer *,complex *, 
	     integer *, complex *,integer *, integer *); 

  float rho2i;
  char s;
  int k,row;
  integer info,n,lda,ldb,nrhs;

  /* Our equation is    y^p_i^* y^q_i r_q = x_i y^p_i*    */
  /* The lapack matrix convention is that the n x n matrix A_pq 
     is stored in a 1-dimensional array of size n^2 with
                A_pq = A[p+n*q]                           
     i.e we cycle fast down the rows of the first column 
     and then slowly across to the next column. Eg packing is as
                 ( 0   3   6 )
		 ( 1   4   7 )
		 ( 2   5   8 )                            */



  /* set lapack parameters for chesv  */

  /*************************************************************/
  /*                                                           */
  /*   Do not change these values without checking that        */
  /*    the changes do not affect the memory allocation        */
  /*                in the calling program                     */
  /*                                                           */
  /*************************************************************/ 


  n=(integer) nenv_chan; /* Size of matrix */
  lda=(integer) nenv_chan;
  ldb=(integer) nenv_chan;
  nrhs=1;    
  s='U';   /* We give the upper triangular part of the matrix */


  
  /*         Call to the lapack routine chesv to solve the equation AR=B        */
  /*   Set up so that R holds the value B on input and the solution on output   */
 
  for (row=0;row<nenv_chan;row++) {
    R[row].r=B[row].r;
    R[row].i=B[row].i;
  }


  /* Call to lapack routine chesv_    */
			                                                  
  chesv_(&s,&n,&nrhs,A,&lda,ipivot,R,&ldb,work,&lwork,&info);
  if (info != 0) {
    GR_start_error("main()",rcsid,__FILE__,__LINE__);
    GR_report_error("Non-zero completion code %d returned by chesv \n",info);
    GR_end_error();  
  }

  /* Calculate the significance measure |rho|^2  */

  *rho2=0.0;
  rho2i=0.0;
  for (row=0;row<nenv_chan;row++) {
    *rho2+=(R[row].r*B[row].r+R[row].i*B[row].i);
    rho2i+=(R[row].i*B[row].r-R[row].r*B[row].i);
  }    
  *rho2/=modx2sum;
  rho2i/=modx2sum;
  if ( fabs(rho2i) > 0.00001 ) {
    GR_start_error("main()",rcsid,__FILE__,__LINE__);
    GR_report_error("Round-off problem:  offset=%d   rho2i=%f (should be 0)\n",offset,rho2i); 
    GR_end_error();  
  }
   
  if ( *rho2 < threshold )  
     for (row=0;row<nenv_chan;row++) {
       R[row].r=R[row].i=0.0;
     }

  /* Calculate the `cleaned' signal  */
    
  for (k=offset;k<offset+correlation_width;k++) {
    rp_clean[k]=rp_signal[k];
    ip_clean[k]=ip_signal[k];
    for (row=0;row<nenv_chan;row++) {
      rp_clean[k]-=(R[row].r*rp_env[row][k]-R[row].i*ip_env[row][k]);
      ip_clean[k]-=(R[row].r*ip_env[row][k]+R[row].i*rp_env[row][k]);
    }
  }
  
  return(0);
}







