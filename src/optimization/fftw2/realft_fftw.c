/* GRASP: Copyright 1997, 1998  Bruce Allen */
/* This is a replacement for the Numerical Recipes realft() routine
   that invokes the FFTW routine fftw.

   Plans are estimated unless the MEASURE flag is set (to do so,
   compile with a -DMEASURE flag).  Note that Wisdom is used with
   measurements---see the FFTW manual for possible problems with this.

   Written by Jolien Creighton <jolien@tapir.caltech.edu>  */

#include "grasp.h"
#include "fftw.h"

static char *rcsid="$Id";

void realft(float data[], unsigned long n, int s)
{
  static int sz = 0;
  static fftw_plan fwd, inv;
  static float *tmp = NULL;
  float *arr = data + 1;
  double theta, wpr, wpi, wr, wi;
  int i;

  if (n != 2*sz) {
    double theta;
    if (sz != 0) {
      fftw_destroy_plan(fwd);
      fftw_destroy_plan(inv);
    }
    sz = n/2;
    tmp = realloc(tmp,n*sizeof(float));
    GR_start_error("realft()",rcsid,__FILE__,__LINE__);
    GR_report_error("using the optimized FFT package FFTW-2\n");
#   ifdef MEASURE
      GR_report_error("measuring plan for n = %d\n",n);
      fwd = fftw_create_plan(sz,FFTW_BACKWARD,FFTW_MEASURE|FFTW_USE_WISDOM);
      inv = fftw_create_plan(sz,FFTW_FORWARD,FFTW_MEASURE|FFTW_USE_WISDOM);
#   else
      GR_report_error("estimating plan for n = %d\n",n);
      fwd = fftw_create_plan(sz,FFTW_BACKWARD,FFTW_ESTIMATE);
      inv = fftw_create_plan(sz,FFTW_FORWARD,FFTW_ESTIMATE);
#   endif
    GR_end_error();
  }

  if (s > 0) {
    memcpy(tmp,arr,n*sizeof(float));
    fftw_one(fwd,(fftw_complex *)tmp,(fftw_complex *)arr);
  }

  theta = (s*M_PI)/(double)sz;
  wpr = cos(theta) - 1;
  wpi = sin(theta);
  wr = 1 + wpr;
  wi = wpi;

  for (i = 1; i < sz/2; ++i) {  /* case i = 0 done separately below */
    int ir = i + i;     /* first half real index */
    int ii = 1 + ir;    /* first half imag index */
    int jr = 2*sz - ir; /* second half real index */
    int ji = jr + 1;    /* second half imag index */
    float h1r =    0.5*(arr[ir] + arr[jr]);
    float h1i =    0.5*(arr[ii] - arr[ji]);
    float h2r =  s*0.5*(arr[ii] + arr[ji]);
    float h2i = -s*0.5*(arr[ir] - arr[jr]);
    arr[ir] =  h1r + wr*h2r - wi*h2i;  /* recombine to form true transform */
    arr[ii] =  h1i + wr*h2i + wi*h2r;
    arr[jr] =  h1r - wr*h2r + wi*h2i;
    arr[ji] = -h1i + wr*h2i + wi*h2r;
    { /* recursion */
      double wt = wr;  /* temporary variable */
      wr = wr*wpr - wi*wpi + wr;  /* recurrence */
      wi = wi*wpr + wt*wpi + wi;
    }
  }

  { /* case i = 0 */
    float h1r = arr[0];
    arr[0] = h1r + arr[1];
    arr[1] = h1r - arr[1];
  }
  if (s < 0) {
    arr[0] *= 0.5;
    arr[1] *= 0.5;
    memcpy(tmp,arr,n*sizeof(float));
    fftw_one(inv,(fftw_complex *)tmp,(fftw_complex *)arr);
  }

  return;
}

