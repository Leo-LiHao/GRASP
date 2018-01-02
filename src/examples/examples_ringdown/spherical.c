/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"

#define TWOPI 6.28318530718
#define FOURPI 12.5663706144
static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))
static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

float sw_spherical(float mu, int s, int l, int m)
/* Computes the spin-weighted spherical harmonic (with phi=0) using
   equation (3.1) of Goldberg et al (1967). */
{
  float factrl(int);
  float bico(int, int);
  float sum,coef,x;
  int sign,r,rmin,rmax;

  if (mu==-1.0) {
    fprintf(stderr,"error in sw_spherical(): division by zero");
    return 0;
  } else {
    x = (1 + mu)/(1 - mu);
  }
  coef = factrl(l+m)*factrl(l-m)*(2*l+1)/(factrl(l-s)*factrl(l+s)*FOURPI);
  rmin = IMAX(0,m-s);
  rmax = IMIN(l-s,l+m);
  sum = 0;
  for (r=rmin;r<=rmax;r++) {
    (((l-r+s)%2)==0) ? (sign = 1) : (sign = -1);
    sum += sign*bico(l-s,r)*bico(l+s,r+s-m)*pow(x,0.5*(2*r+s-m));
  }
  sum *= sqrt(coef)*pow(0.5*(1-mu),l);

  return sum;
}


int main(int argc, char *argv[])
{
  float Sre,Sim,Y,norm=1.0/sqrt(TWOPI),mu=0,dmu=0.02;
  float eigenvalues[4];
  int s,l,m;

  /* process arguments */
  if (argc==4) { /* correct number of arguments */
    s = atoi(argv[1]);
    l = atoi(argv[2]);
    m = atoi(argv[3]);
  } else { /* incorrect number of arguments */
    fprintf(stderr,"usage: spherical s l m\n");
    return 1;
  }

  /* set the eigenvalues to produce spin-weighted spherical harmonics */
  eigenvalues[0] = eigenvalues[1] = eigenvalues[3] = 0;
  eigenvalues[2] = (l - s)*(l + s + 1);

  /* reset the normalization */
  sw_spheroid(&Sre,&Sim,mu,1,0.0,s,l,m,eigenvalues);

  for (mu=-1+0.5*dmu;mu<1;mu+=dmu) {
    /* compute the spin-weighted spheroidal harmonic */
    sw_spheroid(&Sre,&Sim,mu,0,0.0,s,l,m,eigenvalues);
    /* compute the spin-weighted spherical harmonic */
    Y = sw_spherical(mu,s,l,m);
    /* print results with correct normalization for the spheroidal harmonic */
    printf("%e\t%e\t%e\n",mu,norm*Sre,Y);
  }

  return 0;
}
