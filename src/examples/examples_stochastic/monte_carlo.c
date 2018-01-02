/* GRASP: Copyright 1997,1998  Bruce Allen */
/* main program to illustrate monte_carlo() */

#include "grasp.h"
void graphout(float,float,int);

#define DETECTORS_FILE "detectors.dat" /* file containing detector info */
#define SITE1_CHOICE 1       /* identification number for site 1 */
#define SITE2_CHOICE 2       /* identification number for site 2 */
#define FAKE_SB 1            /* 1: simulate stochastic background */
                             /* 0: no stochastic background */
#define FAKE_NOISE1 0        /* 1: simulate detector noise at site 1 */
                             /* 0: no detector noise at site 1 */
#define FAKE_NOISE2 0        /* 1: simulate detector noise at site 2 */
                             /* 0: no detector noise at site 2 */
#define WHITEN_OUT1 1        /* 1: whiten output at site 1 */
                             /* 0: don't whiten output at site 1 */
#define WHITEN_OUT2 1        /* 1: whiten output at site 2 */
                             /* 0: don't whiten output at site 2 */
#define N 65536              /* number of data points */
#define DELTA_T (5.0e-5)     /* sampling period (in sec) */
#define OMEGA_0 (1.0e-3)     /* omega_0 */
#define F_LOW (5.0)          /* minimum frequency (in Hz) */
#define F_HIGH (5.0e3)       /* maximum frequency (in Hz) */
#define NUM_RUNS 5           /* number of runs */

int main(int argc,char **argv)
{
  int i,j,last=0,seed= -17;
  float delta_f,tstart=0.0,time_now;

  float site1_parameters[9],site2_parameters[9];
  char site1_name[100],noise1_file[100],whiten1_file[100];
  char site2_name[100],noise2_file[100],whiten2_file[100];

  double *power1,*power2,*whiten1,*whiten2,*gamma12;
  float *out1,*out2;

  /* ALLOCATE MEMORY */ 
  power1=(double *)malloc((N/2)*sizeof(double));
  power2=(double *)malloc((N/2)*sizeof(double));
  whiten1=(double *)malloc(N*sizeof(double));
  whiten2=(double *)malloc(N*sizeof(double));
  gamma12=(double *)malloc((N/2)*sizeof(double));
  out1=(float *)malloc(N*sizeof(float));
  out2=(float *)malloc(N*sizeof(float));

  /* IDENTITY WHITENING FILTERS (IF WHITEN_OUT1=WHITEN_OUT2=0) */
  for (i=0;i<N/2;i++) {
    whiten1[2*i]=whiten2[2*i]=1.0;
    whiten1[2*i+1]=whiten2[2*i+1]=0.0;
  }

  /* CALL DETECTOR_SITE() TO GET SITE PARAMETER INFORMATION */
  detector_site(DETECTORS_FILE,SITE1_CHOICE,site1_parameters,site1_name,
		noise1_file,whiten1_file);
  detector_site(DETECTORS_FILE,SITE2_CHOICE,site2_parameters,site2_name,
		noise2_file,whiten2_file);

  /* CONSTRUCT NOISE POWER SPECTRA, OVERLAP REDUCTION FUNCTION, AND */
  /* (NON-TRIVIAL) WHITENING FILTERS, IF DESIRED */
  delta_f=(float)(1.0/(N*DELTA_T));
  noise_power(noise1_file,N/2,delta_f,power1);
  noise_power(noise2_file,N/2,delta_f,power2);
  overlap(site1_parameters,site2_parameters,N/2,delta_f,gamma12);
  if (WHITEN_OUT1==1) whiten(whiten1_file,N/2,delta_f,whiten1);
  if (WHITEN_OUT2==1) whiten(whiten2_file,N/2,delta_f,whiten2);

  /* SIMULATE STOCHASTIC BACKGROUND AND/OR DETECTOR NOISE */

  for (j=0;j<NUM_RUNS;j++) {
    monte_carlo(FAKE_SB,FAKE_NOISE1,FAKE_NOISE2,N,DELTA_T,OMEGA_0,F_LOW,F_HIGH,
		gamma12,power1,power2,whiten1,whiten2,out1,out2,&seed);

    /* DISPLAY OUTPUT USING XMGR */
    for (i=0;i<N;i++) {
      time_now=tstart+i*DELTA_T;
      printf("%e\t%e\n",time_now,out1[i]);
    }
    if (j==NUM_RUNS-1) last=1;
    graphout(tstart,tstart+N*DELTA_T,last);

    /* UPDATE TSTART */
    tstart+=N*DELTA_T;

  } /* end for (j=0;j<NUM_RUNS;j++) */

  return 0;
}

void graphout(float xmin,float xmax,int last)
{
  static int first=1;
  printf("&\n");

  if (first) { 

    /* first time we draw plot */
    printf("@doublebuffer true\n"); /* keep display from flashing */
    printf("@focus off\n");          
    printf("@world xmin %e\n",xmin); 
    printf("@world xmax %e\n",xmax); 
    printf("@autoscale yaxes\n");    
    printf("@xaxis label \"t (sec)\"\n"); 
    printf("@title \"Simulated Detector Ouput\"\n"); 
    printf("@subtitle \"(stochastic background--whitened)\"\n");
    printf("@redraw \n");
    if (!last) printf("@kill s0\n"); /* kill set; ready to read again */
    
    first=0;
  }
  else {

    /* other timeOAs we draw plot */
    printf("@world xmin %e\n",xmin); 
    printf("@world xmax %e\n",xmax); 
    printf("@autoscale yaxes\n");    
    if (!last) printf("@kill s0\n"); /* kill set; ready to read again */

  }
 
  return;
}
