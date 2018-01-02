/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"

int main() {
   float tstart,time,srate;
   int i,npoint,code;
   short *data;
   struct fgetinput fgetinput;
   struct fgetoutput fgetoutput;

   /* specify the number of points of output & allocate array */
   npoint=100;
   data=(short *)malloc(sizeof(short)*npoint);
   fgetinput.npoint=npoint;

   /* we only want one channel of data */
   fgetinput.nchan=1;

   /* use the framefiles() function to find it */
   fgetinput.files=framefiles;

   /* allocate space to store channel names */
   fgetinput.chnames=(char **)malloc(fgetinput.nchan*sizeof(char *));

   /* allocate space for data storage location addresses */
   fgetinput.locations=(short **)malloc(fgetinput.nchan*sizeof(short *));

   /* allocate space for numbers of points returned in each channel */
   fgetoutput.npoint=(int *)malloc(fgetinput.nchan*sizeof(int));

   /* allocate space for ratios of channel sample rates */
   fgetoutput.ratios=(int *)malloc(fgetinput.nchan*sizeof(int));

   /* channel name */
   fgetinput.chnames[0]="IFO_DMRO";

   /* set up different cases */
   if (NULL!=getenv("GRASP_REALTIME")) {
      /* don't care if locked */
      fgetinput.inlock=0;
   }
   else {
      /* only locked  */
      fgetinput.inlock=1;
   }

   fgetinput.seek=0;
   fgetinput.calibrate=0;
   fgetinput.locations[0]=data;

   while (1) {
      /* get npoint points of data */
      code=fget_ch(&fgetoutput,&fgetinput);
      tstart=fgetoutput.dt;
      srate=fgetoutput.srate;

      /* if no data remains, exit loop */
      if (code==0) break;
      /* if starting a new locked segment, print banner */
      if (code==1) {
         printf("____________ NEW LOCKED SEGMENT ____________\n\n");
         printf("  Time (sec)\t   IFO output\n");
      }
      /* now output the data */
      for (i=0;i<npoint;i++) {
         time=tstart+i/srate;
         printf("%f\t%d\n",time,(int)data[i]);
      }
   }
   /* close the data files, and return */
   return 0;
}
