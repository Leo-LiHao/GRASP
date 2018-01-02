/* GRASP: Copyright 1997,1998  Bruce Allen */
/* This is a RISKY replacement for the Numerical Recipes realft()
routine, that involkes the IBM ESSL library.  WARNING -- this
replacement routine uses ++array[n] and ++array[n+1] (but first saves
then at the end restores their values) */

#define DEBUG 0

static char *rcsid="$Id: realft_ibm.c,v 1.2 1998/01/23 17:59:43 ballen Exp $\
\n$Name: RELEASE_1_9_8 $";

#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>

/* required so that ESSL routines can return values to specify auxilliary storage sizes */
#define __ESVERR
#include <essl.h>

/* the name of the ESSL error exit routine, to prevent termination on errors */
extern int enotrm();

/* Sets the ESSL error handler to not terminate on errors 2015 or 2030 */
int set_err() {
	int ierno,inoal,inomes,itrace,irange;

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
	errset(&ierno,&inoal,&inomes,&itrace,&enotrm,&irange);

	return 0;
}

void boundary(char *note,int arrayad) {
	if (arrayad%8==0)
		printf("The array %s is aligned on a 8-byte boundary\n",note);
	else if (arrayad%4==0)
		printf("The array %s is aligned on a 4-byte boundary\n",note);
	else if (arrayad%2==0)
		printf("The array %s is aligned on a 2-byte boundary\n",note);
	else
		printf("The array %s is aligned on a 1-byte boundary\n",note);

	return;
}

void realft(float *array,unsigned n,int type) {
	static int loclastnf=0,loclastni=0;
        double auxjunk[16];
	int lastnf,lastni;
	static double *lfaux[2]={NULL,NULL},*liaux[2]={NULL,NULL};
	double *faux[2],*iaux[2];
	static int lfnaux[2]={0,0},linaux[2]={0,0};
	int fnaux[2],inaux[2];
	float saveme[2],scale;
	int i,init,inc2x,inc2y,npoint,numseq,naux3,isign;
	int retcode;

	/* since NR routines expect unit offset arrays */
#if (DEBUG)
	boundary("Un-modified input Array",(int)array);
#endif
	array++;
#if (DEBUG)
	boundary("Input Array",(int)array);
#endif
	if (type!=1) type=-1;

	/* static variables may get overwritten on data segment, so save */
	lastnf=loclastnf;
	lastni=loclastni;
	for (i=0;i<2;i++) faux[i]=lfaux[i];
	for (i=0;i<2;i++) iaux[i]=liaux[i];
	for (i=0;i<2;i++) fnaux[i]=lfnaux[i];
	for (i=0;i<2;i++) inaux[i]=linaux[i];
	saveme[0]=array[n];
	saveme[1]=array[n+1];

	/* create some internal storage for the arrays, if needed */
	if ((lastnf!=n && type==1) || (lastni!=n && type!=1)) {
		fprintf(stderr,"NOTE: Initializing type=%d auxilliary arrays for n=%d\n",type,n);
		fprintf(stderr,"NOTE: Using the RISKY optimized FFT routine realft_ibm.\n");
		fprintf(stderr,"NOTE: This routine may cause segmentation violations or damage the stack!\n");
		fprintf(stderr,"NOTE: If you suspect problems, re-run after linking with a safe version.\n");
		fprintf(stderr,"GRASP version: %s\n",rcsid);

		/* set the error handler to not terminate */
		set_err();

		/* set up parameters to FFT routine */
		init=1;        	/* initialize stored quantities */
		inc2x=0;       	/* stride of input array will be ignored */
       		inc2y=0;       	/* stride of output array will be ignored */
        	npoint=n; 	/* number of points to FFT */
        	numseq=1;      	/* number of sequences to transform */
        	fnaux[0]=15;    /* aux storage required: set lower than required, get value later */
        	fnaux[1]=0;     /* aux storage required: set lower than required, get value later */
        	inaux[0]=14;    /* aux storage required: set lower than required, get value later */
        	inaux[1]=0;     /* aux storage required: set lower than required, get value later */
        	naux3=0;      	/* aux storage required: non-zero only for data not double aligned */

		/* first call returns sizes of auxilliary arrays */
		if (type==1) {
			retcode=srcft(init,array,inc2x,array,inc2y,&npoint,numseq,isign=-1,scale=1.0,
				auxjunk,&fnaux[0],NULL,&fnaux[1],NULL,&naux3);
#if (DEBUG)
			fprintf(stderr,"For FFT type=%d: npoint=%d naux1=%d naux2=%d naux3=%d\n",
				type,npoint,fnaux[0],fnaux[1],naux3);
#endif
		}
		else {
			retcode=scrft(init,array,inc2x,array,inc2y,&npoint,numseq,isign=1,scale=0.5,
				auxjunk,&inaux[0],NULL,&inaux[1],NULL,&naux3);
#if (DEBUG)
			fprintf(stderr,"For FFT type=%d: npoint=%d naux1=%d naux2=%d naux3=%d\n",
				type,npoint,inaux[0],inaux[1],naux3);
#endif
		}
#if (DEBUG)
		fprintf(stderr,"Getting parameters: return code from srcft/scrft is %d\n",retcode);
#endif

		/* check that array size one of those that can be handled by optimized FFT routine. */
		if (npoint!=n) {
			fprintf(stderr,"Optimized FFT unable to handle %d points; need %d instead.  Exiting...\n",
				n,npoint);
			abort();
		}

		/* now allocate or re-allocate auxilliary array storage as needed */
		if (type==1) {
			for (i=0;i<2;i++) faux[i]=(double *)realloc((void *)faux[i],sizeof(double)*fnaux[i]);
			if (faux[0]==NULL || faux[1]==NULL) {
				fprintf(stderr,"Optimized FFT unable to allocate %d + %d doubles.  Exiting...\n",
					fnaux[0],fnaux[1]);
				abort();
			}
#if (DEBUG)
			boundary("Forward auxilliary Array #1",(int)faux[0]);
			boundary("Forward auxilliary Array #2",(int)faux[1]);
#endif
		}
		else {
			for (i=0;i<2;i++) iaux[i]=(double *)realloc(iaux[i],sizeof(double)*inaux[i]);
			if (iaux[0]==NULL || iaux[1]==NULL) {
				fprintf(stderr,"Optimized FFT unable to allocate %d + %d doubles.  Exiting...\n",
					inaux[0],inaux[1]);
				abort();
			}
#if (DEBUG)
			boundary("Inverse auxilliary Array #1",(int)faux[0]);
			boundary("Inverse3 auxilliary Array #2",(int)faux[1]);
#endif
		}

		/* second call initializes the auxilliary storage arrays */
		if (type==1) {
			retcode=srcft(init,NULL,inc2x,NULL,inc2y,&npoint,numseq,isign,scale,
				faux[0],&fnaux[0],faux[1],&fnaux[1],NULL,&naux3);
			lastnf=n;
		}
		else {
			retcode=scrft(init,NULL,inc2x,NULL,inc2y,&npoint,numseq,isign,scale,
				iaux[0],&inaux[0],iaux[1],&inaux[1],NULL,&naux3);
			lastni=n;
		}
#if (DEBUG)
		fprintf(stderr,"Setting auxilliary storage: return code from srcft/scrft is %d\n",retcode);
		fprintf(stderr,"Done with initialization.\n");
#endif
	}

       	init=0;		/* don't initialize */
	inc2x=0;       	/* stride of input array will be ignored */
       	inc2y=0;       	/* stride of output array will be ignored */
        numseq=1;      	/* number of sequences to transform */
        naux3=0;      	/* aux storage required: non-zero only for data not double aligned */

	/* forward transform */
	if (type==1) {
        	isign=-1;       /* sign of exponential for forward transform */
        	scale=1.0;     	/* overall scaling factor */

		/* do FFT of the time-domain data stream */
		retcode=srcft(init,array,inc2x,array,inc2y,&lastnf,numseq,isign,scale,
			faux[0],&fnaux[0],faux[1],&fnaux[1],auxjunk,&naux3);

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

	        isign=1;       /* sign of exponential for inverse transform */
	        scale=0.5;     /* overall scaling factor */

		/* Now do the inverse FFT */
		retcode=scrft(init,array,inc2x,array,inc2y,&lastni,numseq,isign,scale,
			iaux[0],&inaux[0],iaux[1],&inaux[1],auxjunk,&naux3);

	}
#if (DEBUG)
	fprintf(stderr,"Setting auxilliary storage: return code from srcft/scrft is %d\n",retcode);
	fprintf(stderr,"Done with FFT.\n");
#endif   
	/* restore static variables */
	loclastnf=lastnf;
	loclastni=lastni;
	for (i=0;i<2;i++) lfaux[i]=faux[i];
	for (i=0;i<2;i++) liaux[i]=iaux[i];
	for (i=0;i<2;i++) lfnaux[i]=fnaux[i];
	for (i=0;i<2;i++) linaux[i]=inaux[i];
	array[n]=saveme[0];
	array[n+1]=saveme[1];
	return;
}
