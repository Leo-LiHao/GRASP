/* GRASP: Copyright 1997,1998  Bruce Allen */
/* This is a replacement for the Numerical Recipes realft() routine,
   that involkes the FFTW routine rfftw.  Written by Jolien Creighton
   jolien@tapir.caltech.edu  */

static char *rcsid="$Id: realft_fftw.c,v 1.2 1999/07/07 22:08:46 ballen Exp $\
\n$Name: RELEASE_1_9_8 $";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rfftw.h"

/* base-2 logarithm */
unsigned int lb(unsigned long n)
{
  unsigned int p=0;
  while ((n>>p)!=1) p++;
  return p;
}

/* base-2 anti-logarithm = */
unsigned long alb(unsigned int p)
{
  return 1L<<p;
}

void realft(float data[], unsigned long n, int isign)
{
  static rfftw_plan p_for,p_inv;
  static unsigned long lastn=0;
  static float *work=NULL;

  /* since NR routines expect unit offset arrays */
  data++;

  /* (re)allocate memory to working array and (re)create plans */
  if (lastn!=n) {
    fprintf(stderr,"GRASP: realft(): Using the optimized FFT routine realft_fftw-1\n");
    fprintf(stderr,"    which uses the FFTW routine rfftw.\n");
    fprintf(stderr,"    Initializing wave array for n=%d\n",n);
    /* check to see if n is a power of two */
    if (n!=alb(lb(n))) {
      fprintf(stderr,"GRASP: realft(): n = %d must be a power of two\n",n);
      fprintf(stderr,"GRASP version: %s\n",rcsid);
      exit(-1);
    }
    fprintf(stderr,"GRASP version: %s\n",rcsid);
    fflush(stderr);

    /* create plans for forward and inverse transforms */
    p_for = rfftw_create_plan(n,1,FFTW_ESTIMATE,REAL_TO_COMPLEX);
    p_inv = rfftw_create_plan(n,-1,FFTW_ESTIMATE,COMPLEX_TO_REAL);

    /* reallocate memory for working array */
    work = (float *)realloc(work,(size_t)((n+2)*sizeof(float)));
    lastn = n;
  }

  if (isign==1) { /* forward transform */
    rfftw(p_for,1,(FFTW_COMPLEX *)data,1,0,(FFTW_COMPLEX *)work,1,0);
    memcpy((void *)data,(const void *)work,(size_t)(n*sizeof(float)));

    /* put the nyquist entry into the second array element */
    data[1] = work[n];
  } else { /* inverse transform */
    int i;
    memcpy((void *)work,(const void *)data,(size_t)(n*sizeof(float)));

    /* put the nyquist entry into the last array element */
    work[n] = work[1];

    /* zero the imaginary parts of nyquist and DC */
    work[n+1] = work[1] = 0;
    rfftw(p_inv,1,(FFTW_COMPLEX *)work,1,0,(FFTW_COMPLEX *)data,1,0);

    /* multiply by scaling factor to restore NR conventions */
    for (i=0;i<n;i++) data[i] *= 0.5;
  }

  return;
}
