/* GRASP: Copyright 1997,1998  Bruce Allen */
/* This is a RISKY replacement for the Numerical Recipes realft() routine,
that involkes the CLASSPACK library on the paragon.  WARNING -- this
replacement routine uses ++array[n] and ++array[n+1] (but first saves
then at the end restores their values) */

static char *rcsid="$Id: realft_paragon_risky.c,v 1.4 1998/01/23 17:59:45 ballen Exp $\
\n$Name: RELEASE_1_9_8 $";

#include <stdlib.h>
#include <stdio.h>

void realft(float *array,unsigned n,int type) {
	static int loclastn=0;
	static float *locwsave=NULL;
	int i,np2;
	float norm,mnorm;
	float saveme[2];

	/* static variables may get overwritten on data segment, so save */
	int lastn;
	float *wsave;
	wsave=locwsave;
	lastn=loclastn;
	saveme[0]=array[n];
	saveme[1]=array[n+1];

	/* since NR routines expect unit offset arrays */
	array++;

	/* normalization */
	norm=0.5*((float)n);
	mnorm=-1.0*norm;

	/* create some internal storage for the arrays */
	if (lastn!=n) {
		fprintf(stderr,"NOTE: Initializing wave array for n=%d\n",n);
		fprintf(stderr,"NOTE: Using the RISKY optimized FFT routine realft_paragon_risky.\n");
		fprintf(stderr,"NOTE: This routine may cause segmentation violations or damage the stack!\n");
		fprintf(stderr,"NOTE: If you suspect problems, re-run after linking with the safe version.\n");
		fprintf(stderr,"GRASP version: %s\n",rcsid);
		wsave=(float *)realloc(wsave,sizeof(float)*(2*n+4));
		/* initialize wsave (does not touch array) */
		scfft1d(array,n,0,wsave);
		lastn=n;
	}

	/* forward transform */
	if (type==1) {
		/* do FFT of the time-domain data stream */
		scfft1d(array,n,1,wsave);
		/* then re-arrange the data into NR order */
		array[1]=array[n];
		array[n]=0.0;
		for (i=3;i<n;i+=2) array[i]*= -1.0;
	}
	else {
		/* start with data in NR order */
		/* move Nyquist into end location */
		array[n]=array[1];
		/* clear imaginary part of DC, Nyquist */
		array[1]=0.0;
		array[n+1]=0.0;
		/* multiply array odd elements, even elements */
		np2=n+2;
		for (i=1;i<np2;i+=2) array[i]*=mnorm;
		for (i=0;i<np2;i+=2) array[i]*=norm;
		csfft1d(array,n,1,wsave);
	}

	/* restore static variables */
	array[n]=saveme[0];
	array[n+1]=saveme[1];
	locwsave=wsave;
	loclastn=lastn;
	return;
}
		
