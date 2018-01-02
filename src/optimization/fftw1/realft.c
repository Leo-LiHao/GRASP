/* GRASP: Copyright 1997,1998  Bruce Allen */
#include <stdio.h>
#include "rfftw.h"
#define NPOINT  8

void realft(float*,unsigned long,int);
void printvec(char *title,float *f, int n);

main()
{
	float array[NPOINT],real[NPOINT],complex[NPOINT+2];
	rfftw_plan plan_for, plan_inv;
	int i;

	for (i=0;i<NPOINT;i++)
		array[i]=real[i]=i*i-6.0*i+5.0;
	printvec("Initial data",array,NPOINT);
	realft(array-1,NPOINT,1);
	printvec("Forward NR transform",array,NPOINT);

	/* create plans */
	plan_for=rfftw_create_plan(NPOINT,1,FFTW_ESTIMATE,REAL_TO_COMPLEX);
	plan_inv=rfftw_create_plan(NPOINT,-1,FFTW_ESTIMATE,COMPLEX_TO_REAL);
	rfftw(plan_for,1,(FFTW_COMPLEX *)real,1,0,(FFTW_COMPLEX *)complex,1,0);
	printvec("Forward FFTW transform",complex,NPOINT+2);

	rfftw(plan_inv,1,(FFTW_COMPLEX *)complex,1,0,(FFTW_COMPLEX *)real,1,0);
	printvec("Result after doing inverse FFTW transform",real,NPOINT);

	return 0;
}
