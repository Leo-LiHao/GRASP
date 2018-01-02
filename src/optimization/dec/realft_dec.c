/* GRASP: Copyright 1997,1998  Bruce Allen */
/* This is a replacement for the Numerical Recipes realft() routine,
that invokes the DEC Extended Math Library  (dxml). 
This program uses the (preferred) three step process, see
realft_dec_1step.c for an implementation using the one step process */

static char *rcsid="$Id: realft_dec.c,v 1.3 1998/01/23 17:59:40 ballen Exp $\n$Name: RELEASE_1_9_8 $";

#include <stdio.h>
#include <dxmldef.h>

void realft(float *array,unsigned n,int type) {
	static DXML_S_FFT_STRUCTURE work_struct;
        static float *work=NULL;
	static int lastn=0;
	int sz,status,stride;
        int i,nover2;
        float norm;

	if (lastn==0) {
	        fprintf(stderr,"GRASP: realft(): Using the optimized FFT routines sfft_*\n");
	        fprintf(stderr,"       from the DEC Extended Math Library (dxml).\n");
	        fprintf(stderr,"GRASP version: %s\n",rcsid);
		fflush(stderr);  
	}
        /* since NR routines expect unit offset arrays */
        array++;

        /* normalization */
	nover2=n/2;
        norm=(float) nover2;

	/* dxml initializations */
    	sz=n;          /* type recasting */
        status=1;
	stride=1;
	
	if (lastn!=n) {
	        if (lastn!=0) sfft_exit_(&work_struct);
		fprintf(stderr,"Initializing dxml fft data structures for n=%d\n",n);
		fprintf(stderr,"GRASP version: %s\n",rcsid);
		status = sfft_init_(&sz,&work_struct,&stride);
		/* print error message if sfft initialisation fails */ 
		if (status != 0) fprintf(stderr,"dxml sfft initialisation error: error code %d\n\n",status);
		/* create some internal storage for the working array */
                work=(float *)realloc(work,n*sizeof(float));
		lastn=n;
	}

	/* forward transform */
	if (type==1) {
		/* do FFT of the time-domain data stream 
	           output stored in work  */
	  status = sfft_apply_("r","r","f",array,work,&work_struct,&stride); 
	        /* print error message if sfft fails */ 
	  if (status != 0) printf("dxml sfft forward transform error:  error code %d\n\n",status);  
        	/* then re-arrange the data back from dxml 
                   order in work to NR order in array */	
	  for (i=0;i<nover2;i++) array[2*i]=work[i];
          array[1]=work[nover2];
	  for (i=1;i<nover2;i++) array[2*i+1]=(-work[n-i]);
		}
	else {
		/* move data from NR order in array
                   to dxml order in work */
          for (i=0;i<nover2;i++) work[i]=array[2*i];
          work[nover2]=array[1];
	  for (i=1;i<nover2;i++) work[n-i]=(-array[2*i+1]); 
	       /* do FFT of the frequency-domain data stream 
	           output stored in array  */ 
           status = sfft_apply_("r","r","b",work,array,&work_struct,&stride); 
                /* print error message if sfft fails */ 
	  if (status != 0) fprintf(stderr,"dxml sfft backward transform error:  error code %d\n\n",status);             
                /* scale array to agree with NR */
	   for (i=0;i<n;i++) array[i]*=norm; 
	}


        return;
}
		
