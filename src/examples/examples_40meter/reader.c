/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"

int main(){
   FILE *fp;
   short *data;
   float tblock,time,srate;
   int code,num,size=0,count=0,which=100;
   struct ld_binheader bheader;
   struct ld_mainheader mheader;

   /* open the IFO channel for reading */
   fp=grasp_open("GRASP_DATAPATH","channel.0","r");

   /* read the first 20 blocks of lock data */
   while (count <20) {
      /* read a block of data */
      code= read_block(fp,&data,&num,&tblock,&srate,1,&size,0,&bheader,&mheader);

      /* if there is no data left, then break */
      if (code==0) break;

      /* print some information about the data.*/
      printf("Data block %d from file channel.0 starts at t = %f sec.\n",count,tblock);
      printf("This block sampled at %f Hz and contains %d shorts.\n",srate,num);

      /* print out some information about a single data point from block */
      time=tblock+(which-1.0)/srate;
      printf("Data item %d at time %f is %d.\n",which,time,data[which-1]);
      printf("The next block of data contains %d shorts.\n\n",code);

      /* increment count of # of blocks read.*/
      count++;
   }

   /* print information about the largest memory block allocated */
   printf("The largest memory block allocated by read block() was %d shorts long\n",size);

   /* free the array allocated by read_block() */
   free(data);
   return 0;
}
