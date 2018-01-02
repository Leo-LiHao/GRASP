/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"

int main() {
   float tstart,tend,srate,totaltime,begin,end;
   int start_offset,start_block,end_offset,end_block,points,zero=0;
   struct ld_mainheader mh;
   struct ld_binheader bh;
   double doubleutc;
   FILE *fplock;
   time_t unixtime;

   /* Open the file for reading */
   fplock=grasp_open("GRASP_DATAPATH","channel.10","r");

   /* print the absolute start time (in UTC) of the run */
   read_block(fplock,NULL,&zero,&tstart,&srate,0,&zero,1,&bh,&mh);
   doubleutc=mh.epoch_time_sec+0.001*mh.epoch_time_msec;
   printf("Starting time of first frame: %13f Unix-C time\n",doubleutc);
   printf("Starting time of first frame: %13f GPS time\n",doubleutc-UTCTOGPS);
   unixtime=mh.epoch_time_sec;
   printf(" ~    UTC time %s",asctime(utctime(&unixtime)));
   printf(" ~ Unix gmtime %s\n",asctime(gmtime(&unixtime)));

   /* rewind the file pointer */
   rewind(fplock);

   while (1) {

      /* find the next locked section of the data */
      points=find_locked(fplock,&start_offset,
                     &start_block,&end_offset,&end_block,&tstart,&tend,&srate);

      /* if no data remains, then exit */
      if (points==0)
         break;

      /* calculate start and end of lock times */
      begin=tstart+start_offset/srate;
      end=tend+end_offset/srate;
      totaltime=end-begin;

      /* print out info for lock intervals > 30 seconds */ 
      if (totaltime>30.0) {
         printf("Locked from t = %f to %f for %f sec\n",begin,end,totaltime);
         printf("Number of data points is %d\n",points);
         printf("Start block: %d  End block: %d\n",start_block,end_block);
         printf("Start offset: %d End offset %d\n\n",start_offset,end_offset);
      }
   }
   return 0;
}
