/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"

#define IOTA 0.0       /* inclination (radians) */
#define BETA 0.0       /* azimuth (radians) */
#define EPS 0.03       /* fractional mass loss */
#define MASS 50.0      /* mass (solar masses) */
#define SPIN 0.98      /* specific angular momentum */
#define MODE_L 2       /* mode integer l */
#define MODE_M 2       /* mode integer m */
#define SRATE 16000.0  /* sampling rate (Hz) */
#define ATTEN 20.0     /* attenuation leven (dB) */
#define MAX 1024       /* max number of points in waveform */

int main()
{
  float *plus,*cross,t,dt=1/SRATE;
  int i,n;

  /* set arrays to NULL so that memory is allocated in called routines */
  plus = cross = NULL;

  /* generate the waveform function data */
  n = qn_ring(IOTA,BETA,EPS,MASS,SPIN,MODE_L,MODE_M,dt,ATTEN,MAX,&plus,&cross);

  /* output the data */
  for (i=0,t=0;i<n;i++,t+=dt) printf("%e\t%e\t%e\n",t,plus[i],cross[i]);

  return 0;
}
