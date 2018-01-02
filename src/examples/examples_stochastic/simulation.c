/* GRASP: Copyright 1997,1998  Bruce Allen */
/* main program for stochastic background simulation */

#include "grasp.h"

#define DETECTORS_FILE "detectors.dat" /* file containing detector info */  
#define SITE1_CHOICE 1       /* identification number for site 1 */
#define SITE2_CHOICE 2       /* identification number for site 2 */
#define FAKE_SB 1            /* 1: simulate stochastic background */
                             /* 0: stochastic background from real data */
#define FAKE_NOISE1 1        /* 1: simulate detector noise at site 1 */
                             /* 0: detector noise from real data at site 1 */
#define FAKE_NOISE2 1        /* 1: simulate detector noise at site 2 */
                             /* 0: detector noise from real data at site 2 */
#define N 65536              /* number of data points */
#define DELTA_T (5.0e-5)     /* sampling period (in sec) */
#define OMEGA_0 (1.0e-3)     /* omega_0 */
#define F_LOW (0.0)          /* minimum frequency (in Hz) */
#define F_HIGH (1.0e4)       /* maximum frequency (in Hz) */
#define REAL_TIME_NOISE1 0   /* 1: use real-time noise at site 1 */
		             /* 0: use noise information from data file */
#define REAL_TIME_NOISE2 0   /* 1: use real-time noise at site 2 */
		             /* 0: use noise information from data file */
#define NUM_RUNS 2500        /* number of runs (for simulation) */
#define NUM_BINS 200         /* number of bins (for statistics) */

int main(int argc,char **argv)
{
  int    i,pass_test=0,previous_test,runs_completed=0,seed= -17;
  float  delta_f;
  double signal,variance;

  float  site1_parameters[9],site2_parameters[9];
  char   site1_name[100],noise1_file[100],whiten1_file[100];
  char   site2_name[100],noise2_file[100],whiten2_file[100];

  double *generation_power1,*generation_power2;
  double *analysis_power1,*analysis_power2;
  double *whiten1,*whiten2;
  double *gamma12;
  float  *out1,*out2;
  float  *stats;

  /* ALLOCATE MEMORY */ 
  generation_power1=(double *)malloc((N/2)*sizeof(double));
  generation_power2=(double *)malloc((N/2)*sizeof(double));
  analysis_power1=(double *)malloc((N/2)*sizeof(double));
  analysis_power2=(double *)malloc((N/2)*sizeof(double));
  whiten1=(double *)malloc(N*sizeof(double));
  whiten2=(double *)malloc(N*sizeof(double));
  gamma12=(double *)malloc((N/2)*sizeof(double));
  out1=(float *)malloc(N*sizeof(float));
  out2=(float *)malloc(N*sizeof(float));
  stats=(float *)malloc(NUM_RUNS*sizeof(float));

  /* INITIALIZE OUTPUT ARRAYS TO ZERO */
  for (i=0;i<N;i++) out1[i]=out2[i]=0.0;

  /* CALL DETECTOR_SITE() TO GET SITE PARAMETER INFORMATION */
  detector_site(DETECTORS_FILE,SITE1_CHOICE,site1_parameters,site1_name,
		noise1_file,whiten1_file);
  detector_site(DETECTORS_FILE,SITE2_CHOICE,site2_parameters,site2_name,
		noise2_file,whiten2_file);

  /* DISPLAY STOCHASTIC BACKGROUND SIMULATION PARAMETERS */
  printf("\n");
  printf("STOCHASTIC GRAVITATIONAL WAVE BACKGROUND SIMULATION\n");
  printf("\n");
  printf("PARAMETERS:\n");
  printf("Simulated stochastic background (0=no,1=yes): %d\n",FAKE_SB);
  printf("Simulated detector noise at site 1 (0=no,1=yes): %d\n",FAKE_NOISE1);
  printf("Simulated detector noise at site 2 (0=no,1=yes): %d\n",FAKE_NOISE2);
  printf("Real-time noise at site 1 (0=no,1=yes): %d\n", REAL_TIME_NOISE1);
  printf("Real-time noise at site 2 (0=no,1=yes): %d\n", REAL_TIME_NOISE2);
  printf("Detector site 1 = %s\n",site1_name);
  printf("Detector site 2 = %s\n",site2_name);
  printf("Sampling period = %e seconds\n",DELTA_T);
  printf("Number of data points = %d\n",N);
  printf("Omega_0 = %e\n",OMEGA_0);
  printf("f_low  = %e Hz\n",F_LOW);
  printf("f_high = %e Hz\n",F_HIGH);
  printf("Number of runs (for simulation) = %d\n",NUM_RUNS);
  printf("Number of bins (for statistics) = %d\n",NUM_BINS);
  printf("\n");

  /* CONSTRUCT NOISE POWER (FOR SIGNAL GENERATION), WHITENING FILTER */
  /* AND THE OVERLAP REDUCTION FUNCTION */
  delta_f=(float)(1.0/(N*DELTA_T));
  noise_power(noise1_file,N/2,delta_f,generation_power1);
  noise_power(noise2_file,N/2,delta_f,generation_power2);
  whiten(whiten1_file,N/2,delta_f,whiten1);
  whiten(whiten2_file,N/2,delta_f,whiten2);
  overlap(site1_parameters,site2_parameters,N/2,delta_f,gamma12);

  /* CONSTRUCT NOISE_POWER (FOR SIGNAL ANALYSIS) IF REAL-TIME NOISE */
  /* IS NOT DESIRED */
  if (REAL_TIME_NOISE1!=1) {
    for (i=0;i<N/2;i++) analysis_power1[i]=generation_power1[i];
  }
  if (REAL_TIME_NOISE2!=1) {
    for (i=0;i<N/2;i++) analysis_power2[i]=generation_power2[i];
  }

  /* PERFORM THE SIMULATION */
  for (i=1;i<=NUM_RUNS;i++) {

    /* SIMULATE STOCHASTIC BACKGROUND AND/OR DETECTOR NOISE, IF DESIRED */
    if (FAKE_SB==1 || FAKE_NOISE1==1 || FAKE_NOISE2==1) {
      monte_carlo(FAKE_SB,FAKE_NOISE1,FAKE_NOISE2,N,DELTA_T,OMEGA_0,
		  F_LOW,F_HIGH,gamma12,
		  generation_power1,generation_power2,
		  whiten1,whiten2,out1,out2,&seed);
    }

    /* TEST DATA TO SEE IF GAUSSIAN */
    previous_test=pass_test;
    pass_test=test_data12(N,out1,out2);

    if (pass_test==1) {

      /* ANALYZE DATA */
      analyze(previous_test,out1,out2,N,DELTA_T,OMEGA_0,F_LOW,F_HIGH,
	      gamma12,whiten1,whiten2,
	      REAL_TIME_NOISE1,REAL_TIME_NOISE2,
	      analysis_power1,analysis_power2,&signal,&variance);

      /* DISPLAY PRELIMINARY STATISTICS */
      prelim_stats(OMEGA_0,N*DELTA_T,signal,variance);

      /* UPDATE RUNS COMPLETED AND STATS ARRAY FOR FINAL STATISTICS */
      runs_completed++;
      stats[runs_completed-1]=signal;
    }

  } /* end for (i=1;i<+NUM_RUNS;i++) */

  /* FINAL STATISTICS */
  printf("\n");
  statistics(stats,runs_completed,NUM_BINS);

  return 0;
}
