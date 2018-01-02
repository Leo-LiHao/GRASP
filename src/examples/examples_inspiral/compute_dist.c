#include "grasp.h"

#define NPOINTS 32768

int main(int argc,char **argv) {
  void HELP();
  double *S_h;
  float site_parameters[8];
  double delta_f,srate = 20000.; /* Hz */
  char noise_file[128],whiten_file[128],site_name[128];
  int i,npoint=131072;

  /* default parameter values */
  double m1_z=1.4,m2_z=1.4,R=1.e-8,snr=5.5,h100=0.75;
  int detector=21;
  double dinsp,dmerge,zinsp,zmerge,Vcinsp,Vcmerge;

  /* get paramters from the command line */
  for(i=0;i<argc-1;i++) {
    if(!strcmp(argv[i+1],"-h")) {
      HELP();
      exit(0);
    }
    if(!strcmp(argv[i+1],"-snr")) {
      if(!(snr=strtod(argv[i+2],NULL))) {
	fprintf(stderr,"Error assigning SNR, defaulting to 5.5.\n");
	snr=5.5;
      }
      if(snr<0.) snr=5.5;
    }
    if(!strcmp(argv[i+1],"-h100")) {
      if(!(h100=strtod(argv[i+2],NULL))) {
	fprintf(stderr,"Error assigning Hubble constant h100, defaulting to 0.75\n");
	h100=0.75;
      }
      if(h100<0.) h100=fabs(h100);
    }
    if(!strcmp(argv[i+1],"-m1")) {
      if(!(m1_z=strtod(argv[i+2],NULL))) {
	fprintf(stderr,"Error assigning redshifted m1, defaulting to 1.4\n");
	m1_z = 1.4;
      }
      if(m1_z<0.) m1_z=1.4;
    }
    if(!strcmp(argv[i+1],"-m2")) {
      if(!(m2_z=strtod(argv[i+2],NULL))) {
	fprintf(stderr,"Error assigning redshifted m2, defaulting to 1.4\n");
	m2_z = 1.4;
      }
      if(m2_z<0.) m2_z=1.4;
    }
    if(!strcmp(argv[i+1],"-d")) {
      if(!(detector=atoi(argv[i+2]))) {
	fprintf(stderr,"Error assigning detector number, defaulting to 21\n");
	detector=21;
      }
      if(detector<1) detector=21;
    }
    if(!strcmp(argv[i+1],"-R")) {
      if(!(R=strtod(argv[i+2],NULL))) {
	fprintf(stderr,"Error assigning rate density, defaulting to 1.e-8\n");
	R=1.e-8;
      }
      if(R<0.) R=1.e-8;
    }
  }

  /* Get info for that detector */
  detector_site("detectors.dat",detector,site_parameters,site_name,
		noise_file,whiten_file);

  /* allocate memory for the noise power spectrum */
  S_h=(double *)malloc(sizeof(double)*(npoint/2+1));
  delta_f=srate/((double)npoint);

  /* Fill in the noise power spectrum for the detector */
  noise_power(noise_file,npoint/2+1,delta_f,S_h);
  
  /* compute effective distance for which inspiral has SNR = value */
  inspiral_dist(&dinsp,&zinsp,&Vcinsp,m1_z,m2_z,snr,S_h,npoint,srate,h100);

  /* compute effective distance for which merger has SNR = value */
  merger_dist(&dmerge,&zmerge,&Vcmerge,m1_z,m2_z,snr,S_h,npoint,srate,h100);

  printf("%s: D_insp  = %e  z_insp  = %e  Vc_insp  = %e  N = %e\n",
	  site_name,dinsp,zinsp,Vcinsp,R*Vcinsp);

  printf("%s: D_merge = %e  z_merge = %e  Vc_merge = %e  N = %e\n\n",
	  site_name,dmerge,zmerge,Vcmerge,R*Vcmerge);

  return 0;
}

void HELP()
{
  fprintf(stderr,"\nThis GRASP code takes the following flags:\n\n");
  fprintf(stderr," -h:            Show this help message.\n\n");
  fprintf(stderr," -snr [snr]:    Use a signal to noise ratio of snr.\n");
  fprintf(stderr,"                Default value is 5.5.\n\n");
  fprintf(stderr," -m1 [(1+z)m1]: Set redshifted mass 1.\n");
  fprintf(stderr,"                Default value is 1.4.\n\n");
  fprintf(stderr," -m2 [(1+z)m2]: Set redshifted mass 2.\n");
  fprintf(stderr,"                Default value is 1.4.\n\n");
  fprintf(stderr," -d [dn]      : Set detector number to dn.  dn must\n");
  fprintf(stderr,"                be a detector integer defined in the\n");
  fprintf(stderr,"                GRASP file detectors.dat.\n");
  fprintf(stderr,"                Default value is 16, corresponding to\n");
  fprintf(stderr,"                Hanford enhancement 0.\n\n");
  fprintf(stderr," -R [R]       : Set rate density to R.  Units are events\n");
  fprintf(stderr,"                per (year Mpc^3).\n");
  fprintf(stderr,"                Default value is 1.e-8.\n\n");
  fprintf(stderr," -h100 [h100] : Set Hubble constant in units of 100 km/sec/Mpc\n");
  fprintf(stderr,"                Default value is 0.75\n\n");
  return;
}
