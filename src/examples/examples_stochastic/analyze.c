/* GRASP: Copyright 1997,1998  Bruce Allen */
/* example program to demonstrate extract_noise() and extract_signal() */

#include "grasp.h"
void graphout(float,float,int);

#define DETECTORS_FILE "detectors.dat" /* file containing detector info */  
#define SITE1_CHOICE 1       /* identification number for site 1 */
#define SITE2_CHOICE 2       /* identification number for site 2 */
#define N 65536              /* number of data points */
#define DELTA_T (5.0e-5)     /* sampling period (in sec) */
#define OMEGA_0 (1.0e-3)     /* omega_0 */
#define F_LOW (5.0)          /* minimum frequency (in Hz) */
#define F_HIGH (5.0e3)       /* maximum frequency (in Hz) */
#define NUM_RUNS 3           /* number of runs */

main()
{
  int i,j,last=0,seed= -17;
  float delta_f;

  float site1_parameters[9],site2_parameters[9];
  char site1_name[100],noise1_file[100],whiten1_file[100];
  char site2_name[100],noise2_file[100],whiten2_file[100];

  double *power1,*power2,*whiten1,*whiten2,*gamma12;
  float *left1,*left2,*right1,*right2;  
  float *out1,*out2;

  /* ALLOCATE MEMORY */ 
  power1=(double *)malloc((N/2)*sizeof(double));
  power2=(double *)malloc((N/2)*sizeof(double));
  whiten1=(double *)malloc(N*sizeof(double));
  whiten2=(double *)malloc(N*sizeof(double));
  gamma12=(double *)malloc((N/2)*sizeof(double));
  left1=(float *)malloc(N*sizeof(float));
  left2=(float *)malloc(N*sizeof(float));
  right1=(float *)malloc(N*sizeof(float));
  right2=(float *)malloc(N*sizeof(float));
  out1=(float *)malloc(N*sizeof(float));
  out2=(float *)malloc(N*sizeof(float));

  /* CALL DETECTOR_SITE() TO GET SITE PARAMETER INFORMATION */
  detector_site(DETECTORS_FILE,SITE1_CHOICE,site1_parameters,site1_name,
		noise1_file,whiten1_file);
  detector_site(DETECTORS_FILE,SITE2_CHOICE,site2_parameters,site2_name,
		noise2_file,whiten2_file);

  /* CONSTRUCT OVERLAP REDUCTION FUNCTION */
  delta_f=(float)(1.0/(N*DELTA_T));
  overlap(site1_parameters,site2_parameters,N/2,delta_f,gamma12);

  /* identity whitening filter */
  for (i=0;i<N/2;i++) {
    whiten1[2*i]=whiten2[2*i]=1.0;
    whiten1[2*i+1]=whiten2[2*i+1]=0.0;
  }

  /* WHITENING FILTER */
/*
  whiten(whiten1_file,N/2,delta_f,whiten1);
  whiten(whiten2_file,N/2,delta_f,whiten2);
*/
  for (j=0;j<NUM_RUNS;j++) {

    /* SIGNAL GENERATION */
    simulate_sb(N,DELTA_T,OMEGA_0,F_LOW,F_HIGH,gamma12,
		whiten1,whiten2,left1,left2,&seed);
    simulate_sb(N,DELTA_T,OMEGA_0,F_LOW,F_HIGH,gamma12,
		whiten1,whiten2,right1,right2,&seed);

    combine_data(1,N,left1,right1,out1);
    combine_data(2,N,left2,right2,out2);

    /* no combining */
/*
    for (i=0;i<N;i++) out1[i]=left1[i];
*/
    /* SIGNAL ANALYSIS */
    extract_noise(0,1,out1,N,DELTA_T,whiten1,power1);
    extract_noise(0,2,out2,N,DELTA_T,whiten2,power2);

    /* DISPLAY SPECTRA USING XMGR */

    for (i=1;i<N/2;i++)
	if (power1[i]>0.0) printf("%e\t%e\n",i*delta_f,power1[i]);

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
    printf("@yaxis label \"log strain squared (1/Hz)\"\n"); 
    printf("@title \"Auto-Correlation Spectrum\"\n"); 
    printf("@subtitle \"(for a pure stochastic background)\"\n");
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
