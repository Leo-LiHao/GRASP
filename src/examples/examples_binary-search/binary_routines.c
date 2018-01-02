/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"

#define SWAP(a,b) temporary=(*(a));(*(a))=(*(b));(*(b))=temporary;

/* static char *rcsid="$Id: binary_routines.c,v 1.7 1998/06/26 17:11:46 ballen Exp $\n$Name: RELEASE_1_9_8 $"; */

/* Compute the correlation coefficient */
void corr_coef(float *a0, float *a1, float *r, int n,
               float *r00, float *r11, float *r01)
{
  float re0,re1,im0,im1,fac;
  int ncomplex=n/2;

  *r00 = *r11 = *r01 = 0;
  while (ncomplex-->0) {

    re0 = *a0++;
    im0 = *a0++;
    re1 = *a1++;
    im1 = *a1++;
    fac = *r++;

    *r00 += fac*(re0*re0 + im0*im0);
    *r11 += fac*(re1*re1 + im1*im1);
    *r01 += fac*(re0*re1 + im0*im1);

  }
  return;
}

void make_retlifs(float m1, float m2, float *ch1, float *ch2, float fstart,
                  int n, float srate, int *filled, float *t_coal,
                  int err_cd_sprs, int order)
{
  float temporary;
  int i,top,half;
  make_filters(m1,m2,ch1,ch2,fstart,n,srate,filled,t_coal,err_cd_sprs,order);
  top=*filled;
  /* is there an even or odd number of elements? */
  half=(top%2==0)?top/2:(top-1)/2;
  for (i=0;i<half;i++) {
    int j=top-i-1;
    SWAP(ch1+i,ch1+j);
    SWAP(ch2+i,ch2+j);
  }
  return;
}
#undef SWAP
