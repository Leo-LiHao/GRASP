/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"
#define NPOINT 65536

int main() {
   void realft(float*,unsigned long,int);
   float response[NPOINT+2],data[NPOINT],freq;
   float res_real,res_imag,dl_real,dl_imag,c0_real,c0_imag,spectrum,srate,factor;
   int i,npoint;
   short datas[NPOINT];
   struct fgetinput fgetinput;
   struct fgetoutput fgetoutput;

   /* We need only the IFO output */
   fgetinput.nchan=1;

   /* use utility function framefiles() to retrieve file names */
   fgetinput.files=framefiles;

   /* storage for channel names, data locations, points returned, ratios */
   fgetinput.chnames=(char **)malloc(fgetinput.nchan*sizeof(char *));
   fgetinput.locations=(short **)malloc(fgetinput.nchan*sizeof(short *));
   fgetoutput.npoint=(int *)malloc(fgetinput.nchan*sizeof(int));
   fgetoutput.ratios=(int *)malloc(fgetinput.nchan*sizeof(int));

   /* set channel name */
   fgetinput.chnames[0]="IFO_DMRO";


   /* are we in the 40-meter lab? */
   if (NULL!=getenv("GRASP_REALTIME")) {
      /* for Caltech 40-meter lab */
      fgetinput.inlock=0;
   }
   else {
      /* for Nov 1994 data set */
      fgetinput.inlock=1;
   }
   /* number of points to sample and fft (power of 2) */
   fgetinput.npoint=npoint=NPOINT;
   fgetinput.calibrate=1;

   /* the array where we want the data to be put */
   fgetinput.locations[0]=datas;

   /* skip 200 seconds into locked region (just seek, no need for data) */
   fgetinput.seek=1;
   fgetoutput.tstart=fgetoutput.lastlock=0.0;
   while (fgetoutput.tstart-fgetoutput.lastlock<200.0)
      fget_ch(&fgetoutput,&fgetinput);

   /* and get next stretch of data (don't seek, we need data) */
   fgetinput.seek=0;
   fget_ch(&fgetoutput,&fgetinput);

   /* the sample rate */
   srate=fgetoutput.srate;

   /* convert gw signal (ADC counts) from shorts to floats */
   for (i=0;i<NPOINT;i++) data[i]=datas[i];

   /* FFT the data  */
   realft(data-1,npoint,1);

   /* get normalization R(f) using swept sine calibration information from frame */
   GRnormalize(fgetoutput.fri,fgetoutput.frinum,npoint,srate,response);

   /* one-sided power-spectrum normalization, to get meters/rHz */
   factor=sqrt(2.0/(srate*npoint));
   /* compute dl.  Leave off DC (i=0) or Nyquist (i=npoint/2) freq */
   for (i=1;i<npoint/2;i++) {
      /* frequency */
      freq=i*srate/npoint;
      /* real and imaginary parts of tilde c0 */
      c0_real=data[2*i];
      c0_imag=data[2*i+1];
      /* real and imaginary parts of R */
      res_real=response[2*i];
      res_imag=response[2*i+1];
      /* real and imaginary parts of tilde dl */
      dl_real=c0_real*res_real-c0_imag*res_imag;
      dl_imag=c0_real*res_imag+c0_imag*res_real;
      /* |tilde dl| */
      spectrum=factor*sqrt(dl_real*dl_real+dl_imag*dl_imag);
      /* output freq in Hz, noise power in meters/rHz */
      printf("%e\t%e\n",freq,spectrum);
   }
   return 0;
}
