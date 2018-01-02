/* GRASP: Copyright 1997,1998  Bruce Allen */
/* compare_chirps.c: by Benjamin Owen, 1997 */

#include "grasp.h"

#define FSTART 40.      /* GW starting frequency, in Hz                   */
#define SRATE 10000.    /* Sample rate, in Hz                             */
#define LENGTH 32768    /* Twice no of freq bins = samples in time domain */
#define MASS1 10.0      /* mass of first body, in solar masses            */
#define MASS2 10.0      /* mass of second body, in solar masses           */

int main() {
  FILE *fp;
  float t_c, f_c, *sp, *td, *dummy;
  int i,pn_order;
  void realft(float*,unsigned long,int);

  /* Allocate memory for chirps */
  sp = (float *)malloc(sizeof(float)*LENGTH);
  td = (float *)malloc(sizeof(float)*LENGTH);
  dummy = (float *)malloc(sizeof(float)*LENGTH);

  /* order of (v/c) used in chirp calculation comparison */
  pn_order=4;

  /* Generate time-domain chirp for comparison purposes */
  make_filters(MASS1, MASS2, td, dummy, FSTART, LENGTH, SRATE, &i, &t_c, 2000,pn_order);

  /* Generate stationary phase chirp in frequency domain */
  f_c = pow(6,-1.5)/M_PI/(MASS1+MASS2)/TSOLAR;
  sp_filters(MASS1, MASS2, sp, dummy, FSTART, LENGTH, SRATE, f_c, t_c, pn_order);

  /* Transform stationary phase chirp to the time domain */
  realft(sp-1, LENGTH, -1);

  /* Graph both chirps in the time domain.  First output file. */
  fp = fopen("compare_chirps.output", "w");
  for (i=0; i<LENGTH; i++)
    fprintf(fp, "%f %f %f\n", i/SRATE, td[i], sp[i]*2*SRATE/LENGTH);
  fclose(fp);

  /* Now graph the contents of the file using xmgr */
  system("xmgr -nxy compare_chirps.output &");

  return 0;
}
