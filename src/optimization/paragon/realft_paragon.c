/* GRASP: Copyright 1997,1998  Bruce Allen */
/* This is a safe replacement for the Numerical Recipes realft()
routine, that involkes the CLASSPACK library on the paragon. */

static char *rcsid="$Id: realft_paragon.c,v 1.4 1998/01/23 17:59:45 ballen Exp $\
\n$Name: RELEASE_1_9_8 $";

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void realft(float *data,unsigned n,int type) {
	static int lastn=0;
	static float *wsave=NULL;
	static float *array=NULL;
	int i,np2;
	float norm,mnorm;

	/* since NR routines expect unit offset arrays */
	data++;

	/* normalization */
	norm=0.5*((float )n);
	mnorm=-1.0*norm;

	/* create some internal storage for the arrays */
	if (lastn!=n) {
		fprintf(stderr,"Initializing wave array for n=%d\n",n);
		fprintf(stderr,"GRASP version: %s\n",rcsid);
		wsave=(float *)realloc(wsave,sizeof(float)*(2*n+4));
		array=(float *)realloc(array,sizeof(float)*(n+2));
		/* initialize wsave (does not touch array) */
		scfft1d(array,n,0,wsave);
		lastn=n;
	}

	/* forward transform */
	if (type==1) {
		memcpy((void *)array,(const void *)data,(size_t)(sizeof(float)*n));
		/* do FFT of the time-domain data stream */
		scfft1d(array,n,1,wsave);
		/* then re-arrange the data into NR order */
		array[1]=array[n];
		array[n]=0.0;
		for (i=3;i<n;i+=2) array[i]*= -1.0;
		memcpy((void *)data,(const void *)array,(size_t)(sizeof(float)*n));
	}
	else {
		/* start with data in NR order */
		/* move Nyquist into end location */
		memcpy((void *)array,(const void *)data,(size_t)(sizeof(float)*n));
		array[n]=array[1];
		/* clear imaginary part of DC, Nyquist */
		array[1]=0.0;
		array[n+1]=0.0;
		/* multiply array odd elements, even elements */
		np2=n+2;
		for (i=1;i<np2;i+=2) array[i]*=mnorm;
		for (i=0;i<np2;i+=2) array[i]*=norm;
		csfft1d(array,n,1,wsave);
		memcpy((void *)data,(const void *)array,(size_t)(sizeof(float)*n));
	}
	return;
}
		
