/* This is a routine that "drives" get_calibrated_data(), printing out the
   start times of segments as well as the return codes.  It is for testing
   purposes, not for doing anything useful. */

#include "binary_params.h"
#include <stdio.h>
#include <stdlib.h>

/* global variables for passing data from get_calibrated_data() */
double datastart;
float *n_inv_noise,*htilde,*pow_renorm,srate=9868.4208984375;
int num_templates,npoint=NPOINT,new_lock,gauss_test,insert_chirp=INSERT_CHIRP;
int get_calibrated_data(),i;
char filename[128];
FILE *fp;


int main(int argc,char **argv)
{
    int get_calibrated_data(),code,segno=0;
    double tlast=0.0,tstart,tend;
    
    if (argc!=3) 
        {
            fprintf(stderr,"Syntax: %s time1 time2\n where time is sec since Jan 1, 1970\n",argv[0]);
            return 1;
        }
    sscanf(argv[1],"%lf",&tstart);
    sscanf(argv[2],"%lf",&tend);
    printf("Dealing with segments in range of times %f to %f\n",tstart,tend);
    fflush(stdout);
    
/* allocate some memory, assign pointers */
    htilde = (float *)malloc((sizeof(float)*npoint+sizeof(float)*(npoint/2+1)+sizeof(float)));
    n_inv_noise = htilde + npoint;
    pow_renorm = n_inv_noise + npoint/2 + 1;
    
/* now loop, aquiring new segments of data */
    while ((code=get_calibrated_data())) {
        if (datastart>tend) return 0;
        segno++;
        if (datastart >=tstart && datastart<=tend)
            {
                sprintf(filename,"spectrum.%05d",segno);
                fp=fopen(filename,"w");
                fprintf(fp,"# time=%f gauss=%d\n",datastart,gauss_test);
                for (i=0;i<npoint/2+1;i++)
                    fprintf(fp,"%e\n",n_inv_noise[i]);
                fclose(fp);
                
            }
    }
    return 0;    
}
