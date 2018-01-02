/* GRASP: Copyright 1997,1998  Bruce Allen */
/* This is a replacement for the Numerical Recipes realft() routine,
that involkes the SGI computational library. */

static char *rcsid="$Id: realft_sgi.c,v 1.2 1998/01/23 17:59:46 ballen Exp $\n$Name: RELEASE_1_9_8 $";

#include <stdio.h>
#include <fft.h>

void realft(float *array,unsigned n,int type) {
        static int loclastn=0;
        static float *loccoeff=NULL;
        int i,np2;
        float norm;
        float saveme[2];

        /* static variables may get overwritten on data segment, so save */
        int lastn;
        float *coeff;
        coeff=loccoeff;
        lastn=loclastn;
        saveme[0]=array[n];
        saveme[1]=array[n+1];

        /* since NR routines expect unit offset arrays */
        array++;

        /* normalization */
        norm=0.5;

	/* create some internal storage for the arrays */
	if (lastn!=n) {
		fprintf(stderr,"GRASP: realft(): Using the optimized FFT routine scfft\n");
		fprintf(stderr,"    which uses the SGI Computational Library.\n");
		fprintf(stderr,"    Initializing coefficient array for n=%d\n",n);
  		fprintf(stderr,"GRASP version: %s\n",rcsid);
		fflush(stderr);
                coeff=(float *)realloc(coeff,sizeof(float)*(n+15));
		/* initialize coeff (does not touch array) */
		coeff = scfft1dui (n,NULL);
		lastn=n;
	}

	/* forward transform */
	if (type==1) {
		/* do FFT of the time-domain data stream */
                scfft1du( 1, n, array, 1, coeff);
		/* then re-arrange the data into NR order */
		array[1]=array[n];
                array[n]=0.0;
	}
	else {
		/* start with data in NR order */
                /* move Nyquist into end location */
                array[n]=array[1];
                /* clear imaginary part of DC, Nyquist */
                array[1]=0.0;
                array[n+1]=0.0;
                csfft1du ( -1, n, array, 1, coeff);
                /* scale array */
                np2=n+2;
                for (i=0;i<np2;i++) array[i]*=norm;
	}

        /* restore static variables */
        array[n]=saveme[0];
        array[n+1]=saveme[1];
        loccoeff=coeff;
        loclastn=lastn;
        return;
}
		
