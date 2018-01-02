/* This is a routine that "drives" get_calibrated_data(), printing out the
   start times of segments as well as the return codes.  It is for testing
   purposes, not for doing anything useful. */

#include "binary_params.h"
#include <stdio.h>

/* global variables for passing data from get_calibrated_data() */
double datastart;
float *n_inv_noise,*htilde,*pow_renorm,srate=9868.4208984375;
int num_templates,npoint=NPOINT,new_lock,gauss_test,insert_chirp=INSERT_CHIRP;
int get_calibrated_data();

int main()
{
    int get_calibrated_data(),code,segno=0;
    double tlast=0.0;

/* allocate some memory, assign pointers */
    htilde = (float *)malloc((sizeof(float)*npoint+sizeof(float)*(npoint/2+1)+sizeof(float)));
    n_inv_noise = htilde + npoint;
    pow_renorm = n_inv_noise + npoint/2 + 1;
    
/* now loop, aquiring new segments of data */
    while ((code=get_calibrated_data())) {
        if (tlast>datastart) printf("ERROR: times not monotonic!\n");
        tlast=datastart;
        printf("Segment %d start-time %f code %d new_lock %d gauss_test %d\n",segno++,datastart,code,new_lock,gauss_test);
        fflush(stdout);
    }
    return 0;
    
}
