/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"

int main() {
   double begin,end,saveend=0.0;
   struct fgetoutput fgetoutput;
   struct fgetinput fgetinput;
   int firstpass=1;
   time_t unixtime;

   /* this number of samples is about 30 seconds of data */
   fgetinput.npoint=296000;

   /* number of channels needed is one */
   fgetinput.nchan=1;

   /* storage for channel names, data locations, points returned, ratios */
   fgetinput.chnames=(char **)malloc(fgetinput.nchan*sizeof(char *));
   fgetinput.locations=(short **)malloc(fgetinput.nchan*sizeof(short *));
   fgetoutput.npoint=(int *)malloc(fgetinput.nchan*sizeof(int));
   fgetoutput.ratios=(int *)malloc(fgetinput.nchan*sizeof(int));

   /* since we operate in SEEK mode, no space needed for data storage */
   fgetinput.locations[0]=NULL;

   /* use utility function framefiles() to retrieve file names */
   fgetinput.files=framefiles;

   /* get only the locked sections */
   fgetinput.inlock=1;

   /* seek over data (we don't care what the values are!) */
   fgetinput.seek=1;

   /* don't need calibration information */
   fgetinput.calibrate=0;

   /* set channel name  */
   fgetinput.chnames[0]="IFO_DMRO";

   /* start the main loop */
   while (1) {

      /* find the next locked section of the data */
      fget_ch(&fgetoutput,&fgetinput);

      /* print out absolute start time of run */
      if (firstpass) {
         printf("Starting time of first frame: %13f Unix-C time\n",fgetoutput.tfirst);
         printf("Starting time of first frame: %13f GPS time\n",fgetoutput.tfirst_gps);
         unixtime=fgetoutput.tfirst;
         printf("    ~ UTC time %s",asctime(utctime(&unixtime)));
         printf(" ~ Unix gmtime %s\n",asctime(gmtime(&unixtime)));
         firstpass=0;
      }

      /* see if we fell out of lock, and print if we did */
      if (fgetoutput.returnval==1) {

         /* time at whick lock lost (relative to start of run) */
         begin=fgetoutput.lostlock-fgetoutput.tfirst;

         /* time at whick lock aquired (relative to start of run) */
         end=fgetoutput.lastlock-fgetoutput.tfirst;
         if (begin>0.0) {
            printf("In lock from t = %f into run to %f into run for %f sec\n",
                      saveend,begin,begin-saveend);
            printf("Out of lock from t = %f into run to %f into run for %f sec\n",
                      begin,end,end-begin);
         }
         saveend=end;
      }

      /* if no data remains, then exit */
      if (fgetoutput.returnval==0) {
         printf("End of data at time %f\n",fgetoutput.tstart-fgetoutput.tfirst);
         break;
      }
   }
   return 0;
}
