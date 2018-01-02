/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"
#define NPOINT 4096

int main() {
   float cplx[NPOINT+2],srate,freq;
   int npoint,i;
   struct fgetoutput fgetoutput;
   struct fgetinput fgetinput;

   /* we need to ask for some sample values, even though all we want is calibration */
   fgetinput.npoint=256;

   /* number of channels */
   fgetinput.nchan=1;

   /* storage for channel names, data locations, points returned, ratios */
   fgetinput.chnames=(char **)malloc(fgetinput.nchan*sizeof(char *));
   fgetoutput.npoint=(int *)malloc(fgetinput.nchan*sizeof(int));
   fgetoutput.ratios=(int *)malloc(fgetinput.nchan*sizeof(int));

   /* use utility function framefiles() to retrieve file names */
   fgetinput.files=framefiles;

   /* don't care if IFO is in lock */
   fgetinput.inlock=0;

   /* don't need data anyway, so might as well seek */
   fgetinput.seek=1;

   /* but we DO need the calibration information */
   fgetinput.calibrate=1;

   /* set the channel name */
   fgetinput.chnames[0]="IFO_DMRO";

   /* number of points of (imagined) FFT */
   npoint=NPOINT;

   /* now get the data (none) and calibration (what we want) */
   fget_ch(&fgetoutput,&fgetinput);

   /* the fast-channel sample rate */
   srate=fgetoutput.srate;

   /* swept sine calibration array is first argument */
   GRcalibrate(fgetoutput.fri,fgetoutput.frinum,npoint,cplx,srate,2,0);

   /* print out frequency, real, imaginary interpolated values */
   printf("# Freq (Hz)\tReal\t\tImag\n");
   for (i=0;i<=NPOINT/2;i++) {
      freq=i*srate/NPOINT;
      printf("%e\t%e\t%e\n",freq,cplx[2*i],cplx[2*i+1]);
   }
   return 0;
}
