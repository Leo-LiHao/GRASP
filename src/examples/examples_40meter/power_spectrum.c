/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"
#define NPOINT 65536

int main() {
   void realft(float*,unsigned long,int);
   float response[NPOINT+2],data[NPOINT],tstart,freq;
   float res_real,res_imag,dl_real,dl_imag,c0_real,c0_imag,spectrum,srate,factor;
   FILE *fpifo,*fplock,*fpss;
   int i,npoint,remain;
   short datas[NPOINT];

   /* open the IFO output file, lock file, and swept-sine file */
   fpifo=grasp_open("GRASP_DATAPATH","channel.0","r");
   fplock=grasp_open("GRASP_DATAPATH","channel.10","r");
   fpss=grasp_open("GRASP_DATAPATH","swept-sine.ascii","r");

   /* number of points to sample and fft (power of 2) */
   npoint=NPOINT;
   /* skip 200 seconds into locked region (seek=1) */
   while (tstart<200.0)
      get_data(fpifo,fplock,&tstart,npoint,datas,&remain,&srate,1);
   /* and get next stretch of data from TTL locked file (seek=0) */
   get_data(fpifo,fplock,&tstart,npoint,datas,&remain,&srate,0);
   /* convert gw signal (ADC counts) from shorts to floats */
   for (i=0;i<NPOINT;i++) data[i]=datas[i];
   /* FFT the data  */
   realft(data-1,npoint,1);
   /* get normalization R(f) using swept sine file */
   normalize_gw(fpss,npoint,srate,response);
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
