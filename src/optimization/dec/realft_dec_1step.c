/* GRASP: Copyright 1997,1998  Bruce Allen */
/* This is a replacement for the Numerical Recipes realft() routine,
that invokes the DEC Extended Math Library  (dxml). 
This program uses the one step process, see realft_dec.c for an 
implementation using the three step process which is preferred when 
performing many transforms of the same length.*/

static char *rcsid="$Id: realft_dec_1step.c,v 1.3 1998/01/23 17:59:40 ballen Exp $\n$Name: RELEASE_1_9_8 $";

#include <stdio.h>
#include <dxmldef.h>

void realft(float *array,unsigned n,int type) {
        int sz,status,stride;
        int i,nover2;
	float *working_array;
        float norm;

	fprintf(stderr,"GRASP: realft(): Using the optimized FFT routine sfft\n");
	fprintf(stderr,"    which uses the DEC Extended Math Library.\n");
	fprintf(stderr,"GRASP version: %s\n",rcsid);
	fflush(stderr);  

        /* since NR routines expect unit offset arrays */
        array++;

        /* normalization */
	nover2=n/2;
        norm=(float) nover2;

	/* dxml initializations */
        sz=n;
	status=1;
	stride=1;

	/* create some internal storage for the working array */
	working_array=(float *)malloc(n*sizeof(float));

	/* forward transform */
	if (type==1) {
		/* do FFT of the time-domain data stream 
	           output stored in working_array  */
	  status = sfft_("r","r","f",array,working_array,&sz,&stride); 
	        /* print error message if sfft fails */ 
	  if (status != 0) fprintf(stderr,"dxml sfft forward transform error:  error code %d\n\n",
				   status);             
		/* then re-arrange the data back from dxml 
                   order in working_array to NR order in array */	
	  for (i=0;i<nover2;i++) array[2*i]=working_array[i];
          array[1]=working_array[nover2];
	  for (i=1;i<nover2;i++) array[2*i+1]=(-working_array[n-i]);
		}
	else {
		/* move data from NR order in array
                   to dxml order in working_array */
          for (i=0;i<nover2;i++) working_array[i]=array[2*i];
          working_array[nover2]=array[1];
	  for (i=1;i<nover2;i++) working_array[n-i]=(-array[2*i+1]); 
	       /* do FFT of the frequency-domain data stream 
	           output stored in array  */ 
           status = sfft_("r","r","b",working_array,array,&sz,&stride); 
                /* print error message if sfft fails */ 
	   if (status != 0) fprintf(stderr,"dxml sfft backward transform error:  error code %d\n\n",
				    status);             
                /* scale array to agree with NR */
	   for (i=0;i<n;i++) array[i]*=norm; 
	}

	free(working_array);

        return;
}
		
