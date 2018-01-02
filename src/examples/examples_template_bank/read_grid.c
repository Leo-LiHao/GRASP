/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"

int main(int argc, char **argv)
{
  struct cubic_grid grid;
  float m1,m2,coef[10];

  /* Read command line. */
  if(argc!=3){
    fprintf(stderr,"Usage: %s M1 M2\n",argv[0]);
    return 1;
  }
  m1=atof(argv[1]);
  m2=atof(argv[2]);

  /* Read data file. */
  if(read_cubic(&grid,"cubic_coef_40meter_m=0.8-3.2.ascii")){
    fprintf(stderr,"Error: %s: file"
	    " cubic_coef_40meter_m=0.8-3.2.ascii is corrupt.\n",
	    argv[0]);
    return 1;
  }

  /* Get coefficients. */
  if(get_cubic(m1,m2,grid,coef))
    fprintf(stdout,"(%.3f,%.3f) lies outside of grid."
	    "  Extrapolating...\n",m1,m2);

  /* Print coefficients. */
  fprintf(stdout,"Match coefficients: %10.3e %10.3e %10.3e\n",coef[0],
	  coef[1],coef[2]);
  fprintf(stdout,"                    %10.3e %10.3e %10.3e %10.3e\n",
	  coef[3],coef[4],coef[5],coef[6]);
  fprintf(stdout,"Axis lengths:       %10.3e %10.3e\n",coef[7],
	  coef[8]);
  fprintf(stdout,"Axis angle:         %10.3e\n",coef[9]);
  return 0;
}
