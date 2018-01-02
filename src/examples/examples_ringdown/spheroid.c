/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"

#define SPIN 0.98 /* the dimensionless angular momentum parameter */

int main(int argc, char *argv[])
{
  float re,im,mu=0,dmu=0.02,a=SPIN;
  float eigenvalues[4];
  int s,l,m;

  /* process arguments */
  if (argc==4) { /* correct number of arguments */
    s = atoi(argv[1]);
    l = atoi(argv[2]);
    m = atoi(argv[3]);
  } else { /* incorrect number of arguments */
    fprintf(stderr,"usage: spheroid s l m\n");
    return 1;
  }

  /* get the eigenvalues for the appropriate quasinormal mode */
  qn_eigenvalues(eigenvalues,a,s,l,m);

  /* reset the normalization */
  sw_spheroid(&re,&im,mu,1,a,s,l,m,eigenvalues);

  for (mu=-1+0.5*dmu;mu<1;mu+=dmu) {
    /* compute the spin-weighted spheroidal harmonic */
    sw_spheroid(&re,&im,mu,0,0.0,s,l,m,eigenvalues);
    /* print results */
    printf("%e\t%e\t%e\n",mu,re,im);
  }

  return 0;
}
