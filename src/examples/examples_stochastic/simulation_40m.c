/* GRASP: Copyright 1997,1998  Bruce Allen */
/* main program for stochastic background simulation */

#include "grasp.h"

#define DETECTORS_FILE "detectors.dat" 
                             /* file containing detector info */  
#define DATA1_CHANNEL_0 "18nov94.1/channel.0" 
                             /* file containing IFO data (site 1) */
#define DATA2_CHANNEL_0 "19nov94.1/channel.0" 
                             /* file containing IFO data (site 2) */
#define DATA1_CHANNEL_10 "18nov94.1/channel.10" 
                             /* file containing lock info data (site 1) */
#define DATA2_CHANNEL_10 "19nov94.1/channel.10" 
                             /* file containing lock info data (site 2) */
#define SITE1_CHOICE 1       /* identification number for site 1 */
#define SITE2_CHOICE 2       /* identification number for site 2 */
#define FAKE_SB 1            /* 1: simulate stochastic background */
                             /* 0: stochastic background from real data */
#define FAKE_NOISE1 1        /* 1: simulate detector noise at site 1 */
                             /* 0: detector noise from real data at site 1 */
#define FAKE_NOISE2 1        /* 1: simulate detector noise at site 2 */
                             /* 0: detector noise from real data at site 2 */
#define REAL_IFO 0           /* 1: use real interferometer data */
                             /* 0: don't use real interfermeter data */
#define N 65536              /* number of data points */
#define DELTA_T (5.0e-5)     /* sampling period (in sec) */
#define OMEGA_0 (1.0e-3)     /* omega_0 */
#define F_LOW (5.0)          /* minimum frequency (in Hz) */
#define F_HIGH (5.0e3)       /* maximum frequency (in Hz) */
#define REAL_TIME_NOISE1 0   /* 1: use real-time noise at site 1 */
		             /* 0: use noise information from data file */
#define REAL_TIME_NOISE2 0   /* 1: use real-time noise at site 2 */
		             /* 0: use noise information from data file */
#define NUM_RUNS 1600        /* number of runs (for simulation) */
#define NUM_BINS 200         /* number of bins (for statistics) */

/*
   CIT 40 meter parameters 
#define DELTA_T (1.0/9868.420546)
#define OMEGA_0 (300.0)      
#define F_HIGH (2.0e3)       

*/

main()
{
  int i,j,pass_test=0,previous_test,runs_completed=0,seed= -17;
  float delta_f;
  float srate1,srate2;
  double signal,variance;

  float site1_parameters[9],site2_parameters[9];
  char site1_name[100],noise1_file[100],whiten1_file[100];
  char site2_name[100],noise2_file[100],whiten2_file[100];

  double *generation_power1,*generation_power2;
  double *analysis_power1,*analysis_power2;
  double *whiten1,*whiten2;
  double *gamma12;
  float *out_monte_carlo1,*out_monte_carlo2;
  float *out_IFO1,*out_IFO2;
  float *data1,*data2;
  float *stats;

  char   line[100];

  FILE *fp1,*fp2,*fp1lock,*fp2lock;

  /* MEMORY ALLOCATION */ 
  generation_power1=(double *)malloc((N/2)*sizeof(double));
  generation_power2=(double *)malloc((N/2)*sizeof(double));
  analysis_power1=(double *)malloc((N/2)*sizeof(double));
  analysis_power2=(double *)malloc((N/2)*sizeof(double));
  whiten1=(double *)malloc(N*sizeof(double));
  whiten2=(double *)malloc(N*sizeof(double));
  gamma12=(double *)malloc((N/2)*sizeof(double));
  out_monte_carlo1=(float *)malloc(N*sizeof(float));
  out_monte_carlo2=(float *)malloc(N*sizeof(float));
  out_IFO1=(float *)malloc(N*sizeof(float));
  out_IFO2=(float *)malloc(N*sizeof(float));
  data1=(float *)malloc(N*sizeof(float));
  data2=(float *)malloc(N*sizeof(float));
  stats=(float *)malloc(NUM_RUNS*sizeof(float));

  /* INITIALIZE MONTE_CARLO AND IFO OUTPUT ARRAYS TO ZERO */
  for (i=0;i<N;i++) {
    out_monte_carlo1[i]=out_monte_carlo2[i]=0.0;
    out_IFO1[i]=out_IFO2[i]=0.0;
  }

  /* OPEN IFO DATA FILES FOR READING, IF NECESSARY */
  if (REAL_IFO) {
    fp1=grasp_open("GRASP_DATAPATH",DATA1_CHANNEL_0,"r");
    fp2=grasp_open("GRASP_DATAPATH",DATA2_CHANNEL_0,"r");
    fp1lock=grasp_open("GRASP_DATAPATH",DATA1_CHANNEL_10,"r");
    fp2lock=grasp_open("GRASP_DATAPATH",DATA2_CHANNEL_10,"r");
  }

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
  printf("Real interferometer data (0=no,1=yes): %d\n",REAL_IFO);
  printf("Real-time noise at site 1 (0=no,1=yes): %d\n", 
	  REAL_TIME_NOISE1);
  printf("Real-time noise at site 2 (0=no,1=yes): %d\n", 
	  REAL_TIME_NOISE2);
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
/*
  printf("Hit return to start the simulation.\n");
  fgets(line,sizeof(line),stdin);
  printf("\n");
*/

  for (i=1;i<=NUM_RUNS;i++) {

    /* SIMULATE STOCHASTIC BACKGROUND AND/OR DETECTOR NOISE, IF DESIRED */
    if (FAKE_SB==1 || FAKE_NOISE1==1 || FAKE_NOISE2==1) {
      monte_carlo(FAKE_SB,FAKE_NOISE1,FAKE_NOISE2,N,DELTA_T,OMEGA_0,
		  F_LOW,F_HIGH,gamma12,
		  generation_power1,generation_power2,
		  whiten1,whiten2,out_monte_carlo1,out_monte_carlo2,&seed);
    }

    /* USE REAL INTERFEROMETER DATA, IF DESIRED */
    if (REAL_IFO==1) {
      get_IFO12(fp1,fp2,fp1lock,fp2lock2,N,out_IFO1,out_IFO2,
		&srate1,&srate2);
      printf("using real IFO data!\n");
    }

    /* COMBINE SIMULATED AND REAL DATA */
    for (j=0;j<N;j++) {
      data1[j]=out_IFO1[j]+out_monte_carlo1[j];
      data2[j]=out_IFO2[j]+out_monte_carlo2[j];
    }

    /* TEST DATA TO SEE IF GAUSSIAN */
    previous_test=pass_test;
    pass_test=test_data12(N,data1,data2);

    if (pass_test==1) {

      /* ANALYZE DATA */
      analyze(previous_test,data1,data2,N,DELTA_T,OMEGA_0,F_LOW,F_HIGH,
	      gamma12,whiten1,whiten2,
	      REAL_TIME_NOISE1,REAL_TIME_NOISE2,
	      analysis_power1,analysis_power2,
	      &signal,&variance);

      /* PRELIMINARY STATISTICS */
      prelim_stats(OMEGA_0,N*DELTA_T,signal,variance);

      /* update runs_completed and stats array for final statistics */
      runs_completed++;
      stats[runs_completed-1]=signal;
    }

  } /* end for (i=1;i<+NUM_RUNS;i++) */

  /* FINAL STATISTICS */
  printf("\n");
  statistics(stats,runs_completed,NUM_BINS);

  return;
}

