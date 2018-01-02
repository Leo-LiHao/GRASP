/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"
#include <unistd.h>  /* need the header for the sleep() function */

int main() {
   int i,num_points,num_win,num_freq,padded_length,max_lines,num_removed;
   float nwdt,*data,*mtap_spec_init,*mtap_spec_final,freq,fnyquist;
   struct removed_lines *line_list;
   FILE *fpriver;

   /* data length, padded length, num frequencies including DC, Nyquist */
   num_points=395;
   padded_length=1024;
   num_freq=1+padded_length/2;

   /* number of taper windows to use, and time-freq bandwidth */
   num_win=5;
   nwdt=4.0;
   
   /* maximum number of lines to remove */
   max_lines=8;

   /* allocate arrays */
   data= (float *)malloc(sizeof(float)*num_points);
   mtap_spec_init=(float *)malloc(sizeof(float)*num_freq);
   mtap_spec_final=(float *)malloc(sizeof(float)*num_freq);
   line_list=(struct removed_lines *)malloc(sizeof(struct removed_lines)*max_lines);
   
   /* Read Willamette River data from Percival & Walden example, pg 505 */
   fpriver=grasp_open("GRASP_PARAMETERS","willamette_river.dat","r");
   for (i=0;i<395;i++) fscanf(fpriver,"%f",data+i);
   fclose(fpriver);

   /* Since the data is sampled once per month, fnyquist = 6 cyles/year */
   fnyquist=0.5*12;

   /* pop up a graph of the original data */
   graph(data,num_points,1); sleep(5);

   /* now remove the spectral lines from the data set */
   remove_spectral_lines(data,num_points,padded_length,nwdt,num_win,
          max_lines,500,&num_removed,line_list,mtap_spec_init,mtap_spec_final,1,0,num_freq);

   /* pop up a graph of the original multitapered spectrum */
   graph(mtap_spec_init,num_freq,1); sleep(5);

   /* pop up a graph of the line-removed data and multitapered spectrum */
   graph(data,num_points,1); sleep(5);
   graph(mtap_spec_final,num_freq,1); sleep(5);
      
   /* print out a list of lines removed */
   printf("Total number of lines removed: %d\n",num_removed);
   for (i=0;i<num_removed;i++) {
      freq=line_list[i].index*fnyquist/num_freq;
      printf("Removed line of amplitude %f + i %f at freq %f cycles/year\t",
                line_list[i].re,line_list[i].im,freq);
      printf("(F-test value %f)\n",line_list[i].fvalue);
   }
   return 0;
}
