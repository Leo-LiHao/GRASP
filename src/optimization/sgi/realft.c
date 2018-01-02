/* GRASP: Copyright 1997,1998  Bruce Allen */
#include <stdio.h>
#include <fft.h>
#define REAL 0
#define IMAG 1
#define NPOINT 32

void realft(float*,unsigned long,int);
void printvec(char *title,float *f, int n);

main()
{
	float array1[NPOINT+2],array2[NPOINT+2],coeff[NPOINT+15];
	int i;

        for (i=0;i<NPOINT;i++)
                array1[i]=array2[i]=cos(7.0*i)+i*i-6.0*i+5.0;

        printvec("Initial data",array1,NPOINT);

	realft(array2-1,NPOINT,1);
        printvec("Forward NR transform",array2,NPOINT);

        scfft1dui (NPOINT,coeff);
        scfft1du ( 1, NPOINT, array1, 1, coeff);
                array1[1]=array1[NPOINT];
                array1[NPOINT]=NULL;

        printvec("Forward PERFORMANCE transform",array1,NPOINT+2);

                array1[NPOINT]=array1[1];
                array1[1]=NULL;

        csfft1du ( 1, NPOINT, array1, 1, coeff);
        sscal1d (NPOINT, 1./(float)NPOINT, array1, 1);
        printvec("Result after doing inverse PERFORMANCE transform",array1,NPOINT);

	return 0;
}
