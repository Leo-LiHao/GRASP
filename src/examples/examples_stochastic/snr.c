/* GRASP: Copyright 1997,1998  Bruce Allen */
/* main program to calculate the theoretical snr */

#include "grasp.h"

#define DETECTORS_FILE "detectors.dat" /* file containing detector info */
#define SITE1_CHOICE 1        /* 1=LIGO-Hanford site */
#define SITE2_CHOICE 2        /* 2=LIGO-Livingston site */
#define OMEGA_0 (3.0e-6)      /* Omega_0 (for initial detectors) */
#define F_LOW 0.0             /* minimum frequency (in Hz) */
#define F_HIGH (1.0e+4)       /* maximum frequency (in Hz) */
#define T (1.0e+7)            /* total observation time (in sec) */
#define N 40000               /* number of frequency points */
#define DELTA_F 0.25          /* frequency spacing (in Hz) */

int main(int argc,char **argv)
{
  double mean,variance,stddev,snr;

  float  site1_parameters[9],site2_parameters[9];
  char   site1_name[100],noise1_file[100],whiten1_file[100];
  char   site2_name[100],noise2_file[100],whiten2_file[100];

  double *power1,*power2;
  double *gamma12;

  /* ALLOCATE MEMORY */
  power1=(double *)malloc(N*sizeof(double));
  power2=(double *)malloc(N*sizeof(double));
  gamma12=(double *)malloc(N*sizeof(double));

  /* CALL DETECTOR_SITE() TO GET SITE PARAMETER INFORMATION */
  detector_site(DETECTORS_FILE,SITE1_CHOICE,site1_parameters,site1_name,
		noise1_file,whiten1_file);
  detector_site(DETECTORS_FILE,SITE2_CHOICE,site2_parameters,site2_name,
		noise2_file,whiten2_file);

  /* CALL NOISE_POWER() AND OVERLAP() */
  noise_power(noise1_file,N,DELTA_F,power1);
  noise_power(noise2_file,N,DELTA_F,power2);
  overlap(site1_parameters,site2_parameters,N,DELTA_F,gamma12);

  /* CALCULATE MEAN, VARIANCE, STDDEV, AND SNR */
  mean=OMEGA_0*T;
  variance=calculate_var(N,DELTA_F,OMEGA_0,F_LOW,F_HIGH,T,gamma12,
			 power1,power2);
  stddev=sqrt(variance);
  snr=mean/stddev;

  /* DISPLAY RESULTS */
  printf("\n");
  printf("Detector site 1 = %s\n",site1_name);
  printf("Detector site 2 = %s\n",site2_name);
  printf("Omega_0 = %e\n",OMEGA_0);
  printf("f_low  = %e Hz\n",F_LOW);
  printf("f_high = %e Hz\n",F_HIGH);
  printf("Observation time T = %e sec\n",T);
  printf("Theoretical S/N = %e\n",snr);
  printf("\n");

  return 0;
}
