/* GRASP: Copyright 1997,1998  Bruce Allen */
#include <fcntl.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <stdlib.h>
#include <memory.h>

/* number of points in array to use for testing here */
#define NPOINT  16

/* required so that ESSL routines can return values to specify auxilliary storage sizes */
#define __ESVERR
#include <essl.h>

/* the name of the ESSL error exit routine, to prevent termination on errors */
extern int enotrm();

void realft(float*,unsigned long,int);
void printvec(char *title,float *f,int n);

/* Sets the ESSL error handler to not terminate on errors 2015 or 2030 */
int set_err() {
	int ierno,inoal,inomes,itrace,irange,iusadr,dummy;
	einfo(0,&dummy,&dummy);

	/* Error 2015 generated when auxilliary storage lengths not correct */
	ierno=2015;    /* error option number */
	inoal=0;       /* do not change number of errors */
	inomes=-1;     /* suppress all messages */
	itrace=0;      /* ignored, but MUST be 0,1,2 */
	irange=2015;   /* apply these settings only to error code 2015 */
	errset(&ierno,&inoal,&inomes,&itrace,&enotrm,&irange);
		
	/* Error 2030 generated when array storage length not correct */
	ierno=2030;
	inoal=0;
	inomes=-1;
	itrace=0;
	irange=2030;
        iusadr=0;
	errset(&ierno,&inoal,&inomes,&itrace,&enotrm,&irange);
}

main()
{
	float array1[NPOINT],array2[NPOINT+2];
	int i;
	int init,inc2x,inc2y,npoint,numseq,isign,naux1,naux2,naux3;
        double *aux1,*aux2,*aux3;
        float scale;

	/* load an array with some fairly "random" numbers, and print */
	for (i=0;i<NPOINT;i++)
		array1[i]=array2[i]=i*i-6.0*i+5.0;
	printvec("Initial data",array1,NPOINT);

	/* do a forward Numerical Recipes conventions transform, and print */
	realft(array1-1,NPOINT,1);
	printvec("Forward NR transform",array1,NPOINT);

	/* do an inverse Numerical Recipes conventions transform, and print */
	realft(array1-1,NPOINT,-1);
	printvec("Inverse NR transform",array1,NPOINT);

	/* set the error handler to not terminate */
	set_err();

	/* set up parameters to FFT routine */
	init=1;        /* initialize stored quantities */
	inc2x=0;       /* stride of input array will be ignored */
        inc2y=0;       /* stride of output array will be ignored */
        npoint=NPOINT; /* number of points to FFT */
        numseq=1;      /* number of sequences to transform */
        isign=-1;      /* sign of exponential for forward transform */
        scale=1.0;     /* overall scaling factor */
        naux1=15;      /* aux storage required: set lower than required,get value later */
        naux2=0;       /* aux storage required: set lower than required,get value later */
        naux3=0 ;      /* aux storage required: non-zero only for data not double aligned */

	/* first call returns sizes of auxilliary arrays */
	srcft(init,array2,inc2x,array2,inc2y,&npoint,numseq,isign,scale,
		aux1,&naux1,aux2,&naux2,aux3,&naux3);

	printf("Values returned: npoint=%d naux1=%d naux2=%d naux3=%d\n",npoint,naux1,naux2,naux3);
        aux1=(double *)malloc(sizeof(double)*naux1);
        aux2=(double *)malloc(sizeof(double)*naux2);
        aux3=(double *)malloc(sizeof(double)*naux3);

	/* second call initializes the auxilliary storage arrays */
	srcft(init,array2,inc2x,array2,inc2y,&npoint,numseq,isign,scale,
		aux1,&naux1,aux2,&naux2,aux3,&naux3);

	/* third call does the FFT */
        init=0;
	srcft(init,array2,inc2x,array2,inc2y,&npoint,numseq,isign,scale,
		aux1,&naux1,aux2,&naux2,aux3,&naux3);
	printvec("Forward ESSL transform",array2,NPOINT+2);


	/* set up parameters to FFT routine */
	init=1;        /* initialize stored quantities */
	inc2x=0;       /* stride of input array will be ignored */
        inc2y=0;       /* stride of output array will be ignored */
        npoint=NPOINT; /* number of points to FFT */
        numseq=1;      /* number of sequences to transform */
        isign=1;       /* sign of exponential for forward transform */
        scale=0.5;     /* overall scaling factor */
        naux1=14;      /* aux storage required: set lower than required,get value later */
        naux2=0;       /* aux storage required: set lower than required,get value later */
        naux3=0 ;      /* aux storage required: non-zero only for data not double aligned */

	/* first call returns sizes of auxilliary arrays */
	scrft(init,array2,inc2x,array2,inc2y,&npoint,numseq,isign,scale,
		aux1,&naux1,aux2,&naux2,aux3,&naux3);

	printf("Values returned: npoint=%d naux1=%d naux2=%d naux3=%d\n",npoint,naux1,naux2,naux3);
        aux1=(double *)realloc(aux1,sizeof(double)*naux1);
        aux2=(double *)realloc(aux2,sizeof(double)*naux2);
        aux3=(double *)realloc(aux3,sizeof(double)*naux3);


	/* second call initializes the auxilliary storage arrays */
	scrft(init,array2,inc2x,array2,inc2y,&npoint,numseq,isign,scale,
		aux1,&naux1,aux2,&naux2,aux3,&naux3);

	/* Now do the inverse FFT */
        init=0;
	scrft(init,array2,inc2x,array2,inc2y,&npoint,numseq,isign,scale,
		aux1,&naux1,aux2,&naux2,aux3,&naux3);
	printvec("Result after doing inverse ESSL transform",array2,NPOINT);
	return 0;
}
