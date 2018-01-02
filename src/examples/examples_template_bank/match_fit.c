/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"

#define DETECTOR_NUM 15     /* Smooth fit to Caltech 40m prototype */
#define FLO 120.            /* Hz - low frequency cut off for filtering */
#define FTAU 140.           /* Hz - frequency used in definitions of
			            tau0, tau1. */
/*#define DETECTOR_NUM 8       Caltech 40m prototype */
/*#define DETECTOR_NUM 1       LIGO initial interferometer */
/*#define DETECTOR_NUM 12      LIGO Advanced interferometer */

int main(int argc,char **argv)
{
  float *pfit,*cfit,semimajor,semiminor,theta;
  float m1,m2,matchcont;
  float srate=50000;
  float site_parameters[9];
  char noise_file[128],whiten_file[128],site_name[128];
  int order,tstp,tstc;

  /* Check that the program is called with the correct number of
     arguments; print out argument information if it's not. */
  if(argc != 5) {
    fprintf(stderr,"4 Arguments: 1. Mass of body 1 (solar masses)\n");
    fprintf(stderr,"             2. Mass of body 2 (solar masses)\n");
    fprintf(stderr,"             3. Match contour match value;\n");
    fprintf(stderr,"             4. Waveform order [power of (v/c)]\n");
    fprintf(stderr,"\nExample: match_fit 1.2 1.6 .97 4\n");
    exit(0);
  }
  /* Assign arguments to variables */
  m1=atof(argv[1]);
  m2=atof(argv[2]);
  matchcont=atof(argv[3]);
  order=atoi(argv[4]);

  /* Get the file names for the desired noise curve */
  detector_site("detectors.dat",DETECTOR_NUM,site_parameters,
		site_name,noise_file,whiten_file);

  printf("\nEvaluating templates for detector: %s using data from file: \"%s\"\n\n\n",
        site_name,noise_file);

  /* Allocate memory for the coefficients used in the parabolic fits */
  pfit=(float *)malloc(sizeof(float)*3);
  cfit=(float *)malloc(sizeof(float)*7);

  /* Try to find a parabolic fit */
  tstp=match_parab(m1,m2,matchcont,order,srate,FLO,FTAU,noise_file,
		   &semimajor,&semiminor,&theta,pfit);
  if(tstp) {
    printf("Found a parabolic fit to the match around template with\n");
    printf("m1=%f, m2=%f.\n\n",m1,m2);
    printf("Semimajor axis of best fit ellipse:    %e ms\n",
	   semimajor*1.e3);
    printf("Semiminor axis of best fit ellipse:    %e ms\n",
	   semiminor*1.e3);
    printf("Angle between semimajor and tau0 axis: %f rad\n",theta);
    printf("Fit: m = 1 + %e x^2 + %e xy + %e y^2\n",
	   pfit[0],pfit[1],pfit[2]);
    printf("[where x=dtau0, y=dtau1]\n");
  } else
    printf("Unable to find parabolic fit.  Attempting cubic fit.\n");
  
  /* If the parabola failed, try to find a cubic fit */
  if(!tstp) {
    tstc=match_cubic(m1,m2,matchcont,order,srate,FLO,FTAU,noise_file,
		     &semimajor,&semiminor,&theta,cfit);
    if(tstc) {
      printf("Found a cubic fit to the match around template with\n");
      printf("m1=%f, m2=%f.\n\n",m1,m2);
      printf("Using ellipse constructed from parabolic part of cubic.\n");
      printf("Semimajor axis of best fit ellipse:    %e ms\n",
	     semimajor*1.e3);
      printf("Semiminor axis of best fit ellipse:    %e ms\n",
	     semiminor*1.e3);
      printf("Angle between semimajor and tau0 axis: %f rad\n",theta);
      printf("Fit: m = 1 + %e x^2 + %e xy + %e y^2\n",
	     cfit[0],cfit[1],cfit[2]);
      printf("         + %e x^3 + %e y^3 + %e x^2 y + %e xy^2\n",
	     cfit[3],cfit[4],cfit[5],cfit[6]);
      printf("[where x=dtau0, y=dtau1]\n");
    } else {
      printf("Unable to find a cubic fit.  Try looking for a match\n");
      printf("contour at smaller match value; or, give up on this\n");
      printf("mass regime.\n");
    }
  }

  return 0;
}
