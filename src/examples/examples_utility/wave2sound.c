/* GRASP: Copyright 1997,1998,1999  Bruce Allen */
/* Contributed by Scott Hughes hughes@owasco.physics.uiuc.edu */
#include <stdio.h>
#include "grasp.h"
#define NMAX 262144

main(int argc, char **argv)
{
  int i,N;
  double junk, hp[NMAX], hc[NMAX];
  short shp[NMAX], shc[NMAX];
  char hpfile[20], hcfile[20];
  FILE *infile;

  if(argc != 4) {
    fprintf(stderr,"Arguments: 1. Data file\n");
    fprintf(stderr,"           2. .au file basename\n");
    fprintf(stderr,"           3. Number of lines in data file to use.\n");
    exit(0);
  }
  if(!(infile = fopen(argv[1],"r"))) {
    fprintf(stderr,"%s does not exist.  Go away.\n",argv[2]);
    exit(1);
  }
  N = atoi(argv[3]);
  if(N > NMAX) {
    fprintf(stderr,"N = %d is too large, truncating to %d\n",N,NMAX);
    N = NMAX;
  }

  for(i = 0; i < N; i++) {
    fscanf(infile,"%lf %lf %lf",&junk,&hp[i],&hc[i]);
    shp[i] = (short)(1.e6*hp[i]);
    shc[i] = (short)(1.e6*hc[i]);
  }

  sprintf(hpfile,"%s_hp",argv[2]);
  sprintf(hcfile,"%s_hc",argv[2]);
  sound(shp,N,hpfile,0);
  sound(shc,N,hcfile,0);
}
