/* GRASP: Copyright 1997,1998  Bruce Allen */
/* main program for stochastic background simulation */

#include "grasp.h"

#define SITE1_CHOICE 8       /* 1=LIGO-Hanford site */
#define SITE2_CHOICE 14       /* 2=LIGO-Livingston site */
#define N 65536              /* number of data points */
#define DELTA_T (5.0e-5)     /* sampling period (in sec) */
#define DELTA_T (1.0/9868.420546)
#define OMEGA_0 (1.0e-3)     /* omega_0 (for initial detectors) */
#define OMEGA_0 (300.0)     /* omega_0 (for initial detectors) */
#define F_LOW (5.0)          /* minimum frequency (in Hz) */
#define F_HIGH (5.0e3)       /* maximum frequency (in Hz) */
#define F_HIGH (2000.0)       /* maximum frequency (in Hz) */
#define REAL_TIME_NOISE1 0   /* 0=no real-time noise at site 1 */
#define REAL_TIME_NOISE2 0   /* 0=no real-time noise at site 2 */
#define NUM_RUNS 100000       /* number of runs (for simulation) */
#define NUM_BINS 20          /* number of bins (for statistics) */
#define FAKE_NOISE 1         /* 1: use simulated noise, 0: use real data */
#define MAX_DATA   32765     /* magnitude of largest allowed data value */
main()
{
  int    i;
  float  delta_f;
  int    seed= -17;
  float value;
  short *datas;
  int test1,test2,runsdone=0;

  float  site1_parameters[9],site2_parameters[9];
  char   site1_name[100],noise1_file[100],whiten1_file[100];
  char   site2_name[100],noise2_file[100],whiten2_file[100];

  double *generation_power1,*generation_power2;
  double *analysis_power1,*analysis_power2;
  double *whiten1,*whiten2;
  double *gamma12;
  float  *data1,*data2;
  float  *stats;

  char   line[100];

  /* MEMORY ALLOCATION */ 
  generation_power1=(double *)malloc((N/2)*sizeof(double));
  generation_power2=(double *)malloc((N/2)*sizeof(double));
  analysis_power1=(double *)malloc((N/2)*sizeof(double));
  analysis_power2=(double *)malloc((N/2)*sizeof(double));
  whiten1=(double *)malloc(N*sizeof(double));
  whiten2=(double *)malloc(N*sizeof(double));
  gamma12=(double *)malloc((N/2)*sizeof(double));
  data1=(float *)malloc(N*sizeof(float));
  data2=(float *)malloc(N*sizeof(float));
  stats=(float *)malloc(NUM_RUNS*sizeof(float));
  datas=(short *)malloc(N*sizeof(short));

  /* CALL DETECTOR_SITE() TO GET SITE PARAMETER INFORMATION */
  detector_site(SITE1_CHOICE,site1_parameters,site1_name,
		noise1_file,whiten1_file);
  detector_site(SITE2_CHOICE,site2_parameters,site2_name,
		noise2_file,whiten2_file);

  /* DISPLAY STOCHASTIC BACKGROUND SIMULATION PARAMETERS */
  printf("\n");
  printf("STOCHASTIC GRAVITATIONAL WAVE BACKGROUND SIMULATION\n");
  printf("\n");
  printf("PARAMETERS:\n");
  printf("Detector site 1 = %s\n",site1_name);
  printf("Detector site 2 = %s\n",site2_name);
  printf("Sampling period = %e seconds\n",DELTA_T);
  printf("Number of data points = %d\n",N);
  printf("Omega_0 = %e\n",OMEGA_0);
  printf("f_low  = %e Hz\n",F_LOW);
  printf("f_high = %e Hz\n",F_HIGH);
  printf("Real-time noise at site 1 (0=no,1=yes): %d\n", 
	  REAL_TIME_NOISE1);
  printf("Real-time noise at site 2 (0=no,1=yes): %d\n", 
	  REAL_TIME_NOISE2);
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
  printf("Hit return to start the simulation.\n");
  fgets(line,sizeof(line),stdin);
  printf("\n");

  for (i=1;i<=NUM_RUNS;i++) {

    monte_carlo(N,DELTA_T,OMEGA_0,F_LOW,F_HIGH,gamma12,
		generation_power1,generation_power2,
		whiten1,whiten2,data1,data2,&seed,FAKE_NOISE);

    for (i=0;i<N;i++) {
        value=data1[i];
        if (fabs(value)>MAX_DATA) {
	   fprintf(stderr,"Warning: data value too large for 16 bits!\n");
           value=(value>0.0)?MAX_DATA:-MAX_DATA;
        }
        datas[i]=(short)floor(value+0.5);
    }
    test1=is_gaussian(datas,N,-MAX_DATA-1,MAX_DATA+1,1);

    for (i=0;i<N;i++) {
        value=data2[i];
        if (fabs(value)>MAX_DATA) {
	   fprintf(stderr,"Warning: data value too large for 16 bits!\n");
           value=(value>0.0)?MAX_DATA:-MAX_DATA;
        }
        datas[i]=(short)floor(value+0.5);
    }
    test2=is_gaussian(datas,N,-MAX_DATA-1,MAX_DATA+1,1);



    if (test1 && test2) {
	analyze(data1,data2,N,DELTA_T,F_LOW,F_HIGH,gamma12,REAL_TIME_NOISE1,REAL_TIME_NOISE2,
	    analysis_power1,analysis_power2,whiten1,whiten2,stats);
        runsdone++;
    }
    else
        printf("Data segment failed test!  Gaussian?  h1: %d  h2:  %d\n",test1,test2);

  }

  /* STATISTICS */
  printf("\n");
  statistics(stats,runsdone,NUM_BINS);

  return;
}


void testme() {
  float *x;
  int y,z;

  graph(x,y,z);
  return;
}
