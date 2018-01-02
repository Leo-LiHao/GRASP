/* GRASP: Copyright 1997,1998  Bruce Allen */
/* main program to illustrate the function optimal_filter() */

#include "grasp.h"

#define DETECTORS_FILE "detectors.dat" /* file containing detector info */
#define SITE1_CHOICE 1              /* 1=LIGO-Hanford site */
#define SITE2_CHOICE 2              /* 2=LIGO-Livingston site */
#define N 500                       /* number of frequency points */
#define DELTA_F 1.0                 /* frequency spacing (in Hz) */
#define F_LOW 0.0                   /* minimum frequency (in Hz) */
#define F_HIGH (1.0e+4)             /* maximum frequency (in Hz) */
#define OUT_FILE "LIGO_filter.dat"  /* output filename */

int main(int argc,char **argv)
{
  int    i;
  double f;
  double abs_value,max;

  float  site1_parameters[9],site2_parameters[9];
  char   site1_name[100],noise1_file[100],whiten1_file[100];
  char   site2_name[100],noise2_file[100],whiten2_file[100];

  double *power1,*power2;
  double *gamma12;
  double *filter12;

  FILE *fp;
  fp=fopen(OUT_FILE,"w");

  /* ALLOCATE MEMORY */
  power1=(double *)malloc(N*sizeof(double));
  power2=(double *)malloc(N*sizeof(double));
  gamma12=(double *)malloc(N*sizeof(double));
  filter12=(double *)malloc(N*sizeof(double));

  /* CALL DETECTOR_SITE() TO GET SITE PARAMETER INFORMATION */
  detector_site(DETECTORS_FILE,SITE1_CHOICE,site1_parameters,site1_name,
		noise1_file,whiten1_file);
  detector_site(DETECTORS_FILE,SITE2_CHOICE,site2_parameters,site2_name,
		noise2_file,whiten2_file);

  /* CALL NOISE_POWER() AND OVERLAP() */
  noise_power(noise1_file,N,DELTA_F,power1);
  noise_power(noise2_file,N,DELTA_F,power2);
  overlap(site1_parameters,site2_parameters,N,DELTA_F,gamma12);

  /* CALL OPTIMAL_FILTER() AND DETERMINE MAXIMUM ABSOLUTE VALUE */
  optimal_filter(N,DELTA_F,F_LOW,F_HIGH,gamma12,power1,power2,filter12);

  max=0.0;
  for (i=0;i<N;i++) {
    abs_value=fabs(filter12[i]);
    if (abs_value>max) max=abs_value;
  }

  /* WRITE FILTER FUNCTION (NORMALIZED TO 1) TO FILE */
  for (i=0;i<N;i++) {
    f=i*DELTA_F;
    fprintf(fp,"%e %e\n",f,filter12[i]/max);
  }

  fclose(fp);

  return 0;
}
