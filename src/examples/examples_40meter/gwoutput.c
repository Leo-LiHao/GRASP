/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"

int main() {
   float tstart,time,srate;
   int remain,i,npoint,code;
   FILE *fp,*fplock;
   short *data;

   /* open the IFO output file and lock file */
   fp=grasp_open("GRASP_DATAPATH","channel.0","r");
   fplock=grasp_open("GRASP_DATAPATH","channel.10","r");

   /* specify the number of points of output & allocate array */
   npoint=100;
   data=(short *)malloc(sizeof(short)*npoint);

   while (1) {
      /* fill the array with npoint points of data */
      code=get_data(fp,fplock,&tstart,npoint,data,&remain,&srate,0);
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
         printf("%f\t%d\n",time,data[i]);
      }
   }
   /* close the data files, and return */
   fclose(fp);
   fclose(fplock);
   return 0;
}
