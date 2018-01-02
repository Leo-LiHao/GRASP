/* GRASP: Copyright 1997,1998  Bruce Allen */
/* main program to illustrate the function overlap() */

#include "grasp.h"

#define DETECTORS_FILE "detectors.dat" /* file containing detector info */
#define SITE1_CHOICE 1               /* 1=LIGO-Hanford site */
#define SITE2_CHOICE 2               /* 2=LIGO-Livingston site */
#define N 10000                      /* number of frequency points */
#define DELTA_F 1.0                  /* frequency spacing (in Hz) */
#define OUT_FILE "LIGO_overlap.dat"  /* output filename */

int main(int argc,char **argv)
{
  int    i;
  double f;

  float  site1_parameters[9],site2_parameters[9];
  char   site1_name[100],noise1_file[100],whiten1_file[100];
  char   site2_name[100],noise2_file[100],whiten2_file[100];

  double *gamma12;

  FILE *fp;
  fp=fopen(OUT_FILE,"w");

  /* ALLOCATE MEMORY */
  gamma12=(double *)malloc(N*sizeof(double));
  
  /* CALL DETECTOR_SITE() TO GET SITE PARAMETER INFORMATION */
  detector_site(DETECTORS_FILE,SITE1_CHOICE,site1_parameters,site1_name,
		noise1_file,whiten1_file);
  detector_site(DETECTORS_FILE,SITE2_CHOICE,site2_parameters,site2_name,
		noise2_file,whiten2_file);

  /* CALL OVERLAP() AND WRITE DATA TO THE FILE */
  overlap(site1_parameters,site2_parameters,N,DELTA_F,gamma12);

  for (i=0;i<N;i++) {
    f=i*DELTA_F;
    fprintf(fp,"%e %e\n",f,gamma12[i]);
  }

  fclose(fp);

  return 0;
}
