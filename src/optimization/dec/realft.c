/* GRASP: Copyright 1997,1998  Bruce Allen */
#include <stdio.h>
#include <math.h>
#include <dxmldef.h>

#define NPOINT 32

void realft(float*,unsigned long,int);
void printvec(char *title,float *f, int n);

main()
{
	float array1[NPOINT],array2[NPOINT];
	int i,status,sz,stride;

        sz=NPOINT;
	status=1;
	stride=1;

        for (i=0;i<NPOINT;i++)
                array1[i]=array2[i]=cos(7.0*i)+i*i-6.0*i+5.0;

        printvec("Initial data",array1,NPOINT);

	realft(array1-1,NPOINT,1);
        printvec("Forward NR transform",array1,NPOINT);

        realft(array1-1,NPOINT,-1);
        for (i=0;i<NPOINT;i++) array1[i]*=2.0/NPOINT;
        printvec("Inverse NR transform",array1,NPOINT);

        status=sfft_("r","r","f",array2,array2,&sz,&stride);
	printf("forward   sfft_() return code: %d\n", status);
        printvec("Forward PERFORMANCE transform",array2,sz);


        status=sfft_("r","r","b", array2,array2,&sz,&stride);
	printf("inverse    sfft_() return code: %d\n", status); 
        printvec("Result after doing inverse PERFORMANCE transform",array2,sz);

	return 0;
}
