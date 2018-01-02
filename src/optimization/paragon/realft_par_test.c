/* GRASP: Copyright 1997,1998  Bruce Allen */
#include <stdio.h>
#define NPOINT  32

void realft(float *,unsigned long,int);
void printvec(char *title,float *f, int n);

main()
{
	float array1[NPOINT+2],array2[NPOINT+2];
	int i;

	for (i=0;i<NPOINT;i++)
		array1[i]=array2[i]=i*i-6.0*i+5.0;

	printvec("Initial data",array2,NPOINT);

	realft(array2-1,NPOINT,1);
	printvec("Forward realft_paragon transform",array2,NPOINT+2);

 	realft(array2-1,NPOINT,-1);
	printvec("Inverse realft_paragon transform",array2,NPOINT+2);

	return 0;
}
