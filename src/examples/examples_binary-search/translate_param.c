/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"

static char *rcsid="$Id: translate_param.c,v 1.3 1998/02/25 23:44:24 ballen Exp $\n$Name: RELEASE_1_9_8 $";

#define FTAU 140.0 /* Frequency used to define tau0,tau1 (Hz). */
#define THETA 0.0  /* Rotation of output axes (degrees). */

int main(int argc, char **argv)
{
  int i;
  double m1,m2,tau0,tau1,x,y,pf,cost,sint;
  FILE *fpin,*fpout;

  /* Open input and output files.
  char *outfile;
  if(argc!=2){
    fprintf(stderr,"Usage: %s filename \n",argv[0]);
    return 1;
  }
  if((fpin=fopen(argv[1],"r"))==NULL){
    fprintf(stderr,"Error: %s: file %s does not exist.\n %s \n",argv[0],
	    argv[1],rcsid);
    return 1;
  }
  sprintf(outfile,"%s.tau",argv[1]);
  fpout=fopen(outfile,"w"); */

  printf("translate_param version %s\n",rcsid);

  fpin=fopen("templates.ascii","r");
  fpout=fopen("templates.ascii.tau","w");

  /* Read in m1,m2 data, transform to tau0,tau1, and write it. */
  pf=M_PI*FTAU;
  cost=cos(THETA*M_PI/180.0);
  sint=sin(THETA*M_PI/180.0);
  printf("%f %f %f\n",THETA,cost,sint);
  for(i=0;fscanf(fpin,"%lf %lf\n",&m1,&m2)!=EOF;i++){
    tau_of_mass(m1,m2,pf,&tau0,&tau1);
    x=tau0*cost+tau1*sint;
    y=-tau0*sint+tau1*cost;
    fprintf(fpout,"%f %f\n",x,y);
  }
  fprintf(stdout,"%i lines transformed.\n",i);

  /* Close files and exit. */
  fclose(fpin);
  fclose(fpout);
  return 0;
}
