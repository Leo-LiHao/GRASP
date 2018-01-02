/* Demonstrates the audio() and graph_short() with (1) 200 Hz Cosine wave
with a quadratically varying amplitude (2) 5.0-5.0 solar mass inspiral
chirp (3) Lai-Shapiro 1.4 solar mass: Displays graph and plays audio
for each. */

#include "grasp.h"
#include <unistd.h>
#define S_RATE 9600
#define NUM_PTS 10000
int i,chirp_pts;
float *wave,*dummy,wavemax,t_coal;
short *snd;
struct LS_physical_constants phys_const;

int main(){
   /* Allocate arrays */
   wave=(float *)malloc(NUM_PTS*sizeof(float));
   dummy=(float *)malloc(NUM_PTS*sizeof(float));
   snd=(short *)malloc(NUM_PTS*sizeof(float));

   /* Generate a 200 Hz Cosine wave with quadratically varying amplitude */
   for (i=0;i<NUM_PTS;i++) wave[i]=i*(NUM_PTS-i)*cos(2.0*M_PI*200.0*i/9600.0)/NUM_PTS;
   /* Convert wave to shorts, rescaling maximum amplitude to SHRT_MAX-2 */
   wavemax=0;
   for (i=0;i<NUM_PTS;i++) if (fabs(wave[i])>wavemax) wavemax=fabs(wave[i]);
   for (i=0;i<NUM_PTS;i++) snd[i]=(short)((SHRT_MAX-2)*wave[i]/wavemax);
   /* Graph and play waveform, then pause briefly */
   graph_short(snd,NUM_PTS);
   audio(snd,NUM_PTS);
   sleep(2);

   /* Same procedure for an inspiral chirp for a 5.0-5.0 solar mass system. */
   make_filters(5.0,5.0,wave,dummy,64.0,NUM_PTS,S_RATE,&chirp_pts,&t_coal,4000,2);
   wavemax=0;
   for (i=0;i<NUM_PTS;i++) if (fabs(wave[i])>wavemax) wavemax=fabs(wave[i]);
   for (i=0;i<NUM_PTS;i++) snd[i]=(short)((SHRT_MAX-2)*wave[i]/wavemax);
   graph_short(snd,NUM_PTS);
   audio(snd,NUM_PTS);
   sleep(2);

   /* Define parameters for Lai-Shapiro waveform */
   phys_const.mass=1.4;
   phys_const.radius=10.0;
   phys_const.distance=1.0;
   phys_const.fmax=1000.0;
   phys_const.inclination=phys_const.Phi_0=0;
   /* Then generate, display and play the Lai-Shapiro waveform. */
   LS_waveform(wave,phys_const,0.0,0.0,0.125*M_PI,1.0/S_RATE,NUM_PTS);
   wavemax=0;
   for (i=0;i<NUM_PTS;i++) if (fabs(wave[i])>wavemax) wavemax=fabs(wave[i]);
   for (i=0;i<NUM_PTS;i++) snd[i]=(short)((SHRT_MAX-2)*wave[i]/wavemax);
   graph_short(snd,NUM_PTS);
   audio(snd,NUM_PTS);

   return 0;
}
