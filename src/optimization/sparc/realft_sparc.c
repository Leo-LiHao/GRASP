/* GRASP: Copyright 1997,1998  Bruce Allen */
/* This is a replacement for the Numerical Recipes realft() routine,
that involkes the Sun performance library on the SPARC. 
*/

static char *rcsid="$Id: realft_sparc.c,v 1.2 1998/01/23 17:59:47 ballen Exp $\
\n$Name: RELEASE_1_9_8 $";

#include <stdlib.h>
#include <stdio.h>

void realft(float *array,unsigned n,int type) {
	static int lastn=0;
	static int tot=0;
	static float *wsave=NULL;
	int i,np2;
	float norm,mnorm,temp;
	void *memmove(void *s1, const void *s2, size_t n);

	/* since NR routines expect unit offset arrays */
	array++;

	/* normalization */
	norm=0.5*((float)n);
	mnorm=-1.0*norm;

	/* create some internal storage for the arrays */
	if (lastn!=n) {
		fprintf(stderr,"GRASP: realft(): Using the optimized FFT routine realft_sparc\n");
		fprintf(stderr,"    which uses the Sun Performance Library.\n");
		fprintf(stderr,"    Initializing wave array for n=%d\n",n);
		fprintf(stderr,"GRASP version: %s\n",rcsid);
		fflush(stderr);
		wsave=(float *)realloc(wsave,sizeof(float)*(2*n+15));
		/* initialize wsave (does not touch array) */
		rffti(n,wsave);
		lastn=n;
	}

	/* forward transform */
	if (type==1) {
		/* do FFT of the time-domain data stream */
		rfftf(n,array,wsave);
		/* then re-arrange the data into NR order */
		temp=array[n-1];
		memmove(&array[2],&array[1],(n-2)*sizeof(float));
		array[1]=temp;
		for (i=3;i<n;i+=2) array[i]*= -1.0;
	}
	else {
		/* start with data in NR order */
		for (i=3;i<n;i+=2) array[i]*= -0.5;
		temp=0.5*array[1];
		for (i=0;i<n;i+=2) array[i]*=0.5;
		memmove(&array[1],&array[2],(n-2)*sizeof(float));
		array[n-1]=temp;
		rfftb(n,array,wsave);
	}

	return;
}
		
