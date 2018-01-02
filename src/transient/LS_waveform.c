/* GRASP: Copyright 1997,1998  Bruce Allen */
/*
 * This function calculates the waveforms for rapidly rotating nascent
 * neutron stars according to the fitting function found by Lai and
 * Shapiro [ApJ 442, 259 (1995)]. 
 *
 * Authors:  Warren G Anderson,  warren@ricci.phys.uwm.edu
 *           Patrick R Brady,    patrick@tapir.caltech.edu
 *
 *
 */
#include "grasp.h"

static char *rcsid="$Id: LS_waveform.c,v 1.3 1998/06/26 17:14:19 ballen Exp $Name: $";

float Amplitude_LS;

void LS_freq_deriv(float t, float u[], float dudt[]){
   dudt[1]=Amplitude_LS*pow(u[1],5.2)*(u[1]-1);
}

void LS_phas_and_freq(double Phi[], float u[], float fmax, 
		float dt, int n_samples);

/*
 * Define the frequency at which we move from one fitting function to
 * another.  See Lai and Shapiro, 1995.
 */
#define KILOMETER (1.e5)      /* cm   */
#define FMAX_THRESH (415.0)

/* define the function */		 
void LS_waveform(float h[], struct LS_physical_constants phys_const,
           float sky_theta, float sky_phi, float polarization, float dt, 
	   int n_samples)
{
  int i;
  float M14,R10,fmax_bar,A1,Amp,Fplus,Fcross,h0_plus,h0_cross,*u;
  double *Phi;


  /* test that memory has been allocated */
  if (h==NULL) {
          GR_start_error("LS_waveform()",rcsid,__FILE__,__LINE__);
          GR_report_error("First argument can not be NULL pointer!\n");
          GR_end_error();
          abort();
  }
  
  /* define constants used in Lai's notes */
  M14 = phys_const.mass/(1.4);
  R10 = phys_const.radius/10.0;
  fmax_bar = phys_const.fmax*pow(R10,1.5)/sqrt(M14);
  if((fmax_bar <= FMAX_THRESH)) 
    A1=pow((fmax_bar/1756.0),2.7); 
  else
    A1=pow((fmax_bar/1525.0),3);

  /* define frequency and phase vectors */   
  u=(float *)malloc(n_samples*sizeof(float));
  Phi=(double *)malloc(n_samples*sizeof(double));
  
  /* define constants used in calculation */
  Amp = ((phys_const.mass*phys_const.mass)*A1
           /(phys_const.distance*phys_const.radius))
        *((M_SOLAR*M_SOLAR)/(KILOMETER*MPC));
  Amplitude_LS = pow((M_SOLAR)*phys_const.mass/
	      (phys_const.radius*KILOMETER),2.5)
    *(A1*A1/(0.016*0.016))*sqrt(M14)/pow(R10,1.5);

  /* find detector response functions */
  beam_pattern(sky_theta,sky_phi,polarization,&Fplus,&Fcross);

  /* find waveform coefficients */
  h0_plus = Fplus*Amp*(1+ cos(phys_const.inclination)
                           *cos(phys_const.inclination))/2.0;
  h0_cross = Fcross*Amp*cos(phys_const.inclination);

  /* integrate phase and frequency */
  LS_phas_and_freq(Phi,u,phys_const.fmax,dt,n_samples);

  /* calculate waveform */
  for(i=0;i<=n_samples-1;++i)
    h[i]=pow(u[i],2.1)*sqrt(1-u[i])*(h0_plus*cos(Phi[i])
	                                 +h0_cross*sin(Phi[i]));
  
  /* Free memory */
  free(Phi);
  free(u);
}

#undef FMAX_THRESH
