/* GRASP: Copyright 1997,1998  Bruce Allen */
#include <stdio.h>
#define REAL 0
#define IMAG 1
#define NPOINT  8

void realft(float*,unsigned long,int);
void printvec(char *title,float *f, int n);

main()
{
	float array1[NPOINT],array2[NPOINT+2],wsave[2*NPOINT+4];
	int i;

	for (i=0;i<NPOINT;i++)
		array1[i]=array2[i]=i*i-6.0*i+5.0;
	printvec("Initial data",array1,NPOINT);
	realft(array1-1,NPOINT,1);
	printvec("Forward NR transform",array1,NPOINT);

	/* initialize wsave */
	scfft1d(array2,NPOINT,0,wsave);
	scfft1d(array2,NPOINT,1,wsave);
	printvec("Forward CLASSPACK transform",array2,NPOINT+2);

	csfft1d(array2,NPOINT,1,wsave);
	printvec("Result after doing inverse CLASSPACK transform",array2,NPOINT+2);


	return 0;
}
