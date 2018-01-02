/* GRASP: Copyright 1997,1998  Bruce Allen */
#include <stdio.h>
#include <memory.h>
#include "grasp.h"
#define HSCALE 1.e20
#define NBINS 16384

int main() {
   float fstart,srate,tcoal,*c0,*c90,*response;
   int filled,i;
   void realft(float*,unsigned long, int);
   FILE *fp;

   /* allocate memory */
   c0=(float*)malloc(sizeof(float)*NBINS);
   c90=(float*)malloc(sizeof(float)*NBINS);
   response=(float*)malloc(sizeof(float)*(NBINS+1));

   /* set start frequency, sample rate, make chirp */
   make_filters(1.4,1.4,c0,c90,fstart=140.0,NBINS,srate=9868.0,&filled,&tcoal,4000,4);
   printf("Chirp length is %d.\n",filled);

   /* Uncomment this line to see the impulse response of the instrument */
   /* for (i=0;i<NBINS;i++) c0[i]=0.0; c0[100]=1.0; */

   /* put chirps into frequency domain */
   realft(c0-1,NBINS,1);

   /* open file containing calibration data, get response, and scale */
   fp=grasp_open("GRASP_DATAPATH","swept-sine.ascii","r");
   normalize_gw(fp,NBINS,srate,response);
   for (i=0;i<NBINS;i++) response[i]*=HSCALE;

   /* avoid floating point errors in inversion */
   response[0]=response[1]=1.e10;

   /* determine IFO channel0 input which would have produced waveform */
   ratio(c0,c0,response,NBINS/2);

   /* invert FFT */
   realft(c0-1,NBINS,-1);

   /* make a graph showing channel.0 */
   printf("File temp.graph contains channel.0 produced by 2 x 1.4 solar masses.\n");
   graph(c0,NBINS,1);

   return 0;
}





