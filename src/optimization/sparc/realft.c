/* GRASP: Copyright 1997,1998  Bruce Allen */
#include <stdio.h>
#define REAL 0
#define IMAG 1
#define NPOINT  8

void realft(float*,unsigned long,int);
void printvec(char *title,float *f, int n);
void rfftf(long int n,float *r,float *wsave);
void rfftb(long int n,float *r,float *wsave);
void rffti(long int n,float *wsave);


main()
{
	float array1[NPOINT],array2[NPOINT],wsave[2*NPOINT+15];
	int i;

	for (i=0;i<NPOINT;i++)
		array1[i]=array2[i]=i*i-6.0*i+5.0;
	printvec("Initial data",array1,NPOINT);
	realft(array1-1,NPOINT,1);
	printvec("Forward NR transform",array1,NPOINT);

	/* initialize wsave */
	rffti(NPOINT,wsave);
	rfftf(NPOINT,array2,wsave);
	printvec("Forward PERFORMANCE transform",array2,NPOINT);

	rfftb(NPOINT,array2,wsave);
	printvec("Result after doing inverse PERFORMANCE transform",array2,NPOINT);


	return 0;
}
