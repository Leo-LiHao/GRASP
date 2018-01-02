/* GRASP: Copyright 1997,1998  Bruce Allen */
/* main program to calculate the minimum detectable omega_0 */

#include "grasp.h"
 
#define DETECTORS_FILE "detectors.dat" /* file containing detector info */
#define SITE1_CHOICE 1        /* 1=LIGO-Hanford site */
#define SITE2_CHOICE 2        /* 2=LIGO-Livingston site */
#define SNR 1.65              /* 1.65=SNR for 95% confidence */
#define F_LOW 3.0             /* minimum frequency (in Hz) */
#define F_HIGH (1.0e+4)       /* maximum frequency (in Hz) */
#define T (1.0e+7)            /* total observation time (in sec) */
#define N 40000               /* number of frequency points */
#define DELTA_F 0.25          /* frequency spacing (in Hz) */
#define PLOT_WANTED 0

#if PLOT_WANTED
FILE *fp;
#endif


int main(int argc,char **argv)
{
  int    i;
  float  f;

  double factor,f3,f6,f9,f12,p1,p2,g2;
  double int1,int2,int3,int4;
  double a,b,c,omega_0;
  int site1_choice=SITE1_CHOICE,site2_choice=SITE2_CHOICE;

  float  site1_parameters[9],site2_parameters[9];
  char   site1_name[100],noise1_file[100],whiten1_file[100];
  char   site2_name[100],noise2_file[100],whiten2_file[100];

  double *power1,*power2;
  double *gamma12;

  /* ALLOCATE MEMORY */
  power1=(double *)malloc(N*sizeof(double));
  power2=(double *)malloc(N*sizeof(double));
  gamma12=(double *)malloc(N*sizeof(double));

  /* Use detector sites specified on command line */
  if (argc==3) {
	site1_choice=atoi(argv[1]);
	site2_choice=atoi(argv[2]);
  }

  /* CALL DETECTOR_SITE() TO GET SITE PARAMETER INFORMATION */
  detector_site(DETECTORS_FILE,site1_choice,site1_parameters,site1_name,
		noise1_file,whiten1_file);
  detector_site(DETECTORS_FILE,site2_choice,site2_parameters,site2_name,
		noise2_file,whiten2_file);

#if PLOT_WANTED
  /* output file of integrand if needed */
  fp=fopen(site1_name,"w");
#endif

  /* CALL NOISE_POWER() AND OVERLAP() */
  noise_power(noise1_file,N,DELTA_F,power1);
  noise_power(noise2_file,N,DELTA_F,power2);
  overlap(site1_parameters,site2_parameters,N,DELTA_F,gamma12);

  /* CALCULATE INTEGRALS FOR VARIANCE */
  int1=int2=int3=int4=0.0;

  for (i=1;i<N;i++) { /* start sum at i=1 to avoid possible division */
                      /* by 0 (e.g., if f_low=0) */
    f=i*DELTA_F;
    if (F_LOW<=f && f<=F_HIGH) {
      f3=f*f*f;
      f6=f3*f3;
      f9=f6*f3;
      f12=f9*f3;
      g2=gamma12[i]*gamma12[i];
      p1=power1[i];
      p2=power2[i];

#if PLOT_WANTED
fprintf(fp,"%e\t%e\n",f,1.e-80*g2/(f6*p1*p2));
#endif

      int1+=DELTA_F*g2/(f6*p1*p2);
      int2+=DELTA_F*g2/(f9*p1*p1*p2);
      int3+=DELTA_F*g2/(f9*p1*p2*p2);
      int4+=DELTA_F*g2*(1.0+g2)/(f12*p1*p1*p2*p2);
    }
  }

  /* CALCULATE COEFFICIENTS OF QUADRATIC EQUATION */ 
  factor=10.0*M_PI*M_PI/(3.0*HUBBLE*HUBBLE);

  a=(int4/int1-2.0*T*int1/(SNR*SNR))/(factor*factor);
  b=(int2+int3)/(int1*factor);
  c=1.0;

  /* SOLVE THE QUADRATIC */
  omega_0=0.5*(-b-sqrt(b*b-4*a*c))/a;
  
  /* DISPLAY RESULTS */
  printf("\n");
  printf("Detector site 1 = %s\n",site1_name);
  printf("Detector site 2 = %s\n",site2_name);
  printf("S/N ratio = %e\n",SNR);
  printf("f_low  = %e Hz\n",F_LOW);
  printf("f_high = %e Hz\n",F_HIGH);
  printf("Observation time T = %e sec\n",T);
  printf("Minumum Omega_0 h_100^2 = %e\n",omega_0);
  printf("(This corresponds to false alarm rate 5%% and false dismissal rate 50%%.)\n");
  printf("With a 5%% false alarm rate and a 5%% false dismissal rate:\n");
  printf("minumum Omega_0 h_100^2 = %e\n",2.0*omega_0);
  printf("See Allen & Romano, PRD59 (1999) 102001 eqn (4.38) for details.\n");
  printf("\n");

#if PLOT_WANTED
  fclose(fp);
#endif

  return 0;
}
