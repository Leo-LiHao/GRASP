/* GRASP: Copyright 1997,1998  Bruce Allen */

/* 
   Effective distance program.

   This code calculates the effective distance by finding the
   distance to which an which angle averaged SNR is (snr) and
   then multiplies by 1.10 (Kip's fudge factor, with corrections as
   found by Eanna and myself).  It also computes comoving volume for
   the redshift, which is the most appropriate volume element for the
   conversion between event rate density and event rate.

   Scott Hughes, written 4 Sept 1998
                 last modified 25 August 1999
*/

#include "grasp.h"

#define MERGE_EPS 0.1
#define GAMMA_HI 0.13
#define GAMMA_LO 0.02
#define FUDGE 1.10     /* Fudge factor due to Thorne, */
                       /* corrected by Flanagan & Hughes */
                       /* See documentation. */

/* Routine is based on Eq. (2.30) of Flanagan and Hughes; the
   effective distance is the distance at which
   $\langle\rho^2\rangle^{1/2} = 5.5$.

   This function has the inspiral spectrum hard-wired in.  The
   reason is that this allows me to work the redshifted mass into
   it in a very simple way. */

void inspiral_dist(double *deff, double *z, double *Vc, double m1_z,
		   double m2_z, double snr, double S_h[], int npoint,
		   double srate,double h100)
{
  double prefactor;
  double integral=0.0;
  double mu_z,M_z;
  double df,f,flo;
  double z_of_D(double D,double h100);
  double Vc_of_z(double z,double h100);
  register int i;

  M_z = TSOLAR*(m1_z + m2_z);
  mu_z = TSOLAR*(m1_z*m2_z/(m1_z+m2_z));
  prefactor = (2./15.)*mu_z*pow(M_z,2./3.)/pow(M_PI,4./3.);

  /* (i=0) == DC --- skip it, since C and my grade school teachers
     don't particularly like dividing by zero. 

     Stop the integral at the f = 0.02/M cutoff. */
 
  flo = GAMMA_LO/M_z;
  f=0.;
  i=1;
  df = srate/npoint;
  while(f < flo && i <= npoint/2) {
    f = ((double)i)*df;
    integral += df*pow(f,-7./3.)/S_h[i];
    i++;
  }
  integral *= prefactor;

  *deff = (FUDGE*sqrt(integral)/snr)*(C_LIGHT/MPC);
  *z = z_of_D(*deff,h100);
  *Vc = Vc_of_z(*z,h100);
}
  
/* Routine is based on Eq. (2.30) of Flanagan and Hughes Phys Rev 57 pg
   4535 (1998); the effective distance is the distance at which
   $\langle\rho^2\rangle^{1/2} = 5.5$.

   This function has the merger spectrum hard-wired in.  The reason is
   that this allows me to work the redshifted mass into it in a very
   simple way. */

void merger_dist(double *deff, double *z, double *Vc, double m1_z,
		 double m2_z, double snr, double S_h[], int npoint,
		 double srate,double h100)
{
  double prefactor,flanhughes_F,integral=0.0;
  double mu_z,M_z;
  double df,f,flo,fhi;
  double z_of_D(double D,double h100);
  double Vc_of_z(double z,double h100);
  register int i;

  M_z = TSOLAR*(m1_z + m2_z);
  mu_z = TSOLAR*(m1_z*m2_z/(m1_z+m2_z));

  flo = GAMMA_LO/M_z;
  fhi = GAMMA_HI/M_z;

  flanhughes_F = (4.*mu_z/M_z)*(4.*mu_z/M_z);

  prefactor = (2./5.)*flanhughes_F*MERGE_EPS*M_z*M_z/
    (M_PI*M_PI*(GAMMA_HI - GAMMA_LO));

  /* (i=0) == DC --- skip it, since C and my grade school teachers
     don't particularly like dividing by zero. */
 
  f = flo;
  i = (int)(npoint*flo/srate);
  df = srate/npoint;
  while(f < fhi && i <= npoint/2) {
    f = ((double)i)*df;
    if (f >= flo)
      integral += df/(f*f*S_h[i]);
    i++;
  }
  integral *= prefactor;

  *deff = (FUDGE*sqrt(integral)/snr)*(C_LIGHT/MPC);
  *z = z_of_D(*deff,h100);
  *Vc = Vc_of_z(*z,h100);

  return;
}

/*
   Formula is from Draza's paper [PRD {\bf 48} 4738 (1993)].
   Assumes the input luminosity distance is in Mpc.
*/
double z_of_D(double D,double h100)
{
  double Omega = 1.;        /* Cosmology choice, no doubt wrong ... */
  double z;

  D *= MPC/C_LIGHT;         /* Now it's in seconds */
  z = -1. + 0.5*Omega + 0.5*HUBBLE*h100*D*Omega
    + 0.5*sqrt((2.*HUBBLE*h100*D+1.)*(Omega-2.)*(Omega-2.));
  return z;
}

/*
  Comoving volume as a function of redshift.  Formula is
  from unpublished notes "Distance measures in cosmology"
  by David W. Hogg.  I have specialized to Omega = 1 in
  this formula.
*/
double Vc_of_z(double z,double h100)
{
  double Vc;
  double sqrt1pz;
  double num,den;

  sqrt1pz = sqrt(1. + z);
  num = -3.*z*z + 14.*(-1. + sqrt1pz) + z*(-16. + 9.*sqrt1pz);
  den = pow(1. + z,1.5)*(-1. + sqrt1pz)*(-1. + sqrt1pz);
  Vc = 1. + num/den;
  Vc *= 32.*M_PI/(3.*HUBBLE*HUBBLE*HUBBLE*h100*h100*h100);
  Vc *= (C_LIGHT/MPC)*(C_LIGHT/MPC)*(C_LIGHT/MPC);
  return(Vc);
}
