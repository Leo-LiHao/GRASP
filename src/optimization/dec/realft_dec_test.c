/* GRASP: Copyright 1997,1998  Bruce Allen */
#include <stdio.h>
#include <math.h>  

#define NPOINT 32 

void realft(float *,unsigned long,int); 
void printvec(char *title,float *f, int n);

main()
{
	float array[NPOINT];
	int i;


	for (i=0;i<NPOINT;i++)
		array[i]=cos(7.0*i)+i*i-6.0*i+5.0;

	printvec("Initial data",array,NPOINT); 

	realft(array-1,NPOINT,1); 
	printvec("Forward realft() transform",array,NPOINT);

 	realft(array-1,NPOINT,-1); 
	printvec("Inverse realft() transform",array,NPOINT); 

	return 0;
}
