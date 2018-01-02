/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"
#define NPOINT 4096

int main() {
	float cplx[NPOINT+2],srate,freq;
	int npoint,i;
	FILE *fp;

	/* open the swept-sine calibration file */
        fp=grasp_open("GRASP_DATAPATH","swept-sine.ascii","r");

	/* number of points of (imagined) FFT */
	npoint=NPOINT;

	/* a sample rate often used for fast channels */
	srate=9868.4208984375;

	/* swept sine calibration filename is first argument */
	calibrate(fp,npoint,cplx,srate,2,0);

	/* print out frequency, real, imaginary interpolated values */
	printf("# Freq (Hz)\tReal\t\tImag\n");
	for (i=0;i<=NPOINT/2;i++) {
		freq=i*srate/NPOINT;
		printf("%e\t%e\t%e\n",freq,cplx[2*i],cplx[2*i+1]);
	}
	return 0;
}
