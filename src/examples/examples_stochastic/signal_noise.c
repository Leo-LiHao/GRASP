/* GRASP: Copyright 1997,1998  Bruce Allen */
/* main program to illustrate extract_noise() and extract_signal() */

#include "grasp.h"
void graphout(float,float,int);

#define DETECTORS_FILE "detectors.dat" /* file containing detector info */  
#define SITE1_CHOICE 1       /* identification number for site 1 */
#define SITE2_CHOICE 2       /* identification number for site 2 */
#define FAKE_SB 0            /* 1: simulate stochastic background */
                             /* 0: no stochastic background */
#define FAKE_NOISE1 1        /* 1: simulate detector noise at site 1 */
                             /* 0: no detector noise at site 1 */
#define FAKE_NOISE2 1        /* 1: simulate detector noise at site 2 */
                             /* 0: no detector noise at site 2 */
#define N 65536              /* number of data points */
#define DELTA_T (5.0e-5)     /* sampling period (in sec) */
#define OMEGA_0 (1.0e-3)     /* omega_0 */
#define F_LOW (5.0)          /* minimum frequency (in Hz) */
#define F_HIGH (5.0e3)       /* maximum frequency (in Hz) */
#define REAL_TIME_NOISE1 1   /* 1: use real-time noise at site 1 */
		             /* 0: use noise information from data file */
#define REAL_TIME_NOISE2 1   /* 1: use real-time noise at site 2 */
		             /* 0: use noise information from data file */
#define NUM_RUNS 5           /* number of runs */

main()
{
  int i,j,last=0,seed= -17;
  float f,delta_f;

  float site1_parameters[9],site2_parameters[9];
  char site1_name[100],noise1_file[100],whiten1_file[100];
  char site2_name[100],noise2_file[100],whiten2_file[100];

  double *generation_power1,*generation_power2;
  double *analysis_power1,*analysis_power2;
  double *whiten1,*whiten2;
  double *gamma12,*signal12;
  float *out1,*out2;

  /* ALLOCATE MEMORY */ 
  generation_power1=(double *)malloc((N/2)*sizeof(double));
  generation_power2=(double *)malloc((N/2)*sizeof(double));
  analysis_power1=(double *)malloc((N/2)*sizeof(double));
  analysis_power2=(double *)malloc((N/2)*sizeof(double));
  whiten1=(double *)malloc(N*sizeof(double));
  whiten2=(double *)malloc(N*sizeof(double));
  gamma12=(double *)malloc((N/2)*sizeof(double));
  signal12=(double *)malloc((N/2)*sizeof(double));
  out1=(float *)malloc(N*sizeof(float));
  out2=(float *)malloc(N*sizeof(float));

  /* CALL DETECTOR_SITE() TO GET SITE PARAMETER INFORMATION */
  detector_site(DETECTORS_FILE,SITE1_CHOICE,site1_parameters,site1_name,
		noise1_file,whiten1_file);
  detector_site(DETECTORS_FILE,SITE2_CHOICE,site2_parameters,site2_name,
		noise2_file,whiten2_file);

  /* CONSTRUCT NOISE POWER (FOR SIGNAL GENERATION), WHITENING */
  /* FILTERS, AND OVERLAP REDUCTION FUNCTION */
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
  
  /* SIMULATE STOCHASTIC BACKGROUND AND/OR DETECTOR NOISE */

  for (j=0;j<NUM_RUNS;j++) {

    /* SIGNAL GENERATION */
    monte_carlo(FAKE_SB,FAKE_NOISE1,FAKE_NOISE2,N,DELTA_T,OMEGA_0,
		F_LOW,F_HIGH,gamma12,
		generation_power1,generation_power2,
		whiten1,whiten2,out1,out2,&seed);

    /* SIGNAL ANALYSIS */
    /* (noise power spectra and cross-correlation spectrum) */
    if (REAL_TIME_NOISE1==1)
      extract_noise(1,1,out1,N,DELTA_T,whiten1,analysis_power1);
    if (REAL_TIME_NOISE2==1)
      extract_noise(1,2,out2,N,DELTA_T,whiten2,analysis_power2);

    extract_signal(1,out1,out2,N,DELTA_T,whiten1,whiten2,signal12);

    /* DISPLAY SPECTRA USING XMGR */

    for (i=1;i<N/2;i++) {
      f=i*delta_f;
      printf("%e\t%e\n",f,analysis_power1[i]);

      /* NOTE: uncomment the following command to display the */
      /* cross-correlation spectrum.  but make sure you also */
      /* change logxy to logx in graphout() to allow for y < 0 */
      /* printf("%e\t%e\n",f,signal12[i]); */ 

    }
    if (j==NUM_RUNS-1) last=1;
    graphout(1.0,delta_f*N/2,last);

  } /* end for (j=0;j<NUM_RUNS;j++) */

  return;
}
void graphout(float xmin,float xmax,int last)
{
  static int first=1;
  printf("&\n");

  if (first) { 

    /* first time we draw plot */
    printf("@doublebuffer true\n"); /* keep display from flashing */
    printf("@focus off\n");          
    printf("@g0 type logxy\n");
    printf("@world xmin %e\n",xmin); 
    printf("@world xmax %e\n",xmax); 
    printf("@xaxis tick minor 1\n");
    printf("@xaxis tick major 2\n");
    printf("@autoscale yaxes\n");    
    printf("@xaxis label \"log f (Hz)\"\n"); 
    printf("@yaxis label \"log P(f) (1/Hz)\"\n"); 
    printf("@title \"Noise Power Spectrum\"\n"); 
    printf("@subtitle \"(for an initial LIGO detector)\"\n");
    printf("@redraw \n");
    if (!last) printf("@kill s0\n"); /* kill set; ready to read again */
    
    first=0;
  }
  else {

    /* other times we draw plot */
    printf("@g0 type logxy\n");
    printf("@world xmin %e\n",xmin); 
    printf("@world xmax %e\n",xmax); 
    printf("@xaxis tick minor 1\n");
    printf("@xaxis tick major 2\n");
    printf("@autoscale yaxes\n");    
    if (!last) printf("@kill s0\n"); /* kill set; ready to read again */

  }
 
  return;
}
