/* GRASP: Copyright 1997,1998  Bruce Allen */
/* main program to illustrate the function overlap() for all detector sites */

#include "grasp.h"

#define DETECTORS_FILE "detectors.dat" /* file containing detector info */
#define N 1000                         /* number of frequency points */
#define DELTA_F 1.0                    /* frequency spacing (in Hz) */
#define NUM_SITES 14                   /* output filename */

main()
{
  int    i,j,k;
  double f;

  float  site1_parameters[9],site2_parameters[9];
  char   site1_name[100],noise1_file[100],whiten1_file[100];
  char   site2_name[100],noise2_file[100],whiten2_file[100];
  char   out_file[100];

  double *gamma12;

  FILE *fp;

  /* MEMORY ALLOCATION */
  gamma12=(double *)malloc(N*sizeof(double));

  for (i=1;i<=NUM_SITES;i++) {
    for (j=1;j<i;j++) {

      /* CALL DETECTOR_SITE() TO GET SITE PARAMETER INFORMATION */
      detector_site(DETECTORS_FILE,i,site1_parameters,site1_name,
		    noise1_file,whiten1_file);
      detector_site(DETECTORS_FILE,j,site2_parameters,site2_name,
		    noise2_file,whiten2_file);

      /* CALL OVERLAP() AND WRITE DATA TO THE FILE */
      overlap(site1_parameters,site2_parameters,N,DELTA_F,gamma12);

      /* WRITE DATA TO A FILE LABELED BY SITE ID NUMBERS */
      sprintf(out_file,"overlap%03d%03d.dat",j,i);
      fp=fopen(out_file,"w");

      /* WRITE HEADER INFO FOR XMGR */
      fprintf(fp,"@world xmin %f\n",0.0);
      fprintf(fp,"@world xmax %f\n",log10((double)N*DELTA_F));
      fprintf(fp,"@autoscale yaxes\n");
      fprintf(fp,"@xaxis label \"Log Freq (Hz) \"\n"); 
      fprintf(fp,"@title \"Overlap Reduction Function\"\n");
      fprintf(fp,"@subtitle  \"(for %s and %s detector sites)\"\n",
	      site1_name,site2_name);
	
      for (k=1;k<N;k++) {
	f=k*DELTA_F;
	fprintf(fp,"%e %e\n",log10(f),gamma12[k]);
      }
      fclose(fp);

    }
  }

  return;
}
