/* GRASP: Copyright 1997,1998  Bruce Allen */
/*
 * This program is an example of how one uses the LS_waveform
 * routine. It simply sets initial values, allocates memory for the wave
 * form, and then allows LS_waveform to do the work!
 */

#include "grasp.h"


int main()
{
  int i,n_samps;
  float *h,dt;
  struct LS_physical_constants phys_const;
  float sky_theta,sky_phi,polarization;

  /* Length and sample rate for filter waveform. */  
  n_samps = 1048576;
  dt = 0.0001;

  /* Beam pattern parameters. */  
  sky_theta = 0.0;
  sky_phi = 0.0;
  polarization = 0.125*M_PI;

  /* Physical parameters of system whose waveform is being modelled. */
  phys_const.mass = 1.4; /* solar masses */
  phys_const.radius = 10.0; /* kilometers */
  phys_const.distance = 10.0; /* megaparsecs */
  phys_const.fmax = 800.0; /* Hz */
  phys_const.inclination = 0.0;
  phys_const.Phi_0 = 0.0;

  /* Allocate waveform memory. */
  h=(float *)malloc(n_samps*sizeof(float));

  /* Produce waveform */
  LS_waveform(h,phys_const,sky_theta,sky_phi,polarization,dt,n_samps);

  /* Output waveform. */
  for(i=0;i<n_samps;i++)
  printf("%f\t%e\n",i*dt,h[i]); 
 
return 0;
}
