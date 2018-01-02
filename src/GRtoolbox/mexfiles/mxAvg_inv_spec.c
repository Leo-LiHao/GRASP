/* compile: mex -g -I/LIGO/GRASP/include mxAvg_inv_spec.c -L/LIGO/GRASP/lib -lgrasp -L/usr/local/lib -lrecipes_c	*/
/*															*/	
/* example run: .m		 							                                */
/*															*/
/* Steve Drasco 													*/
/* Summer 1998														*/

#include "mex.h"
#include "grasp.h"
 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double	decay, *norm;
        float   flo, srate, *htilde, *mean_pow_spec, *twice_inv_noise;
        int     n, i;

        /* mex checks */
        if (nrhs != 6)
                mexErrMsgTxt("Need 6 inputs (flo, srate, n, decay, norm, htilde)");
        if (nlhs != 3)
                mexErrMsgTxt("Need 3 outputs (mean_pow_spec, twice_inv_noise, norm)");

        /* make GRASP input */
	flo = (float) mxGetScalar(prhs[0]);
	srate = (float) mxGetScalar(prhs[1]);
	n = (int) mxGetScalar(prhs[2]);
	decay = (double) mxGetScalar(prhs[3]);
	*norm = (double) mxGetScalar(prhs[4]);
	htilde = (float *) mxMalloc(n * sizeof(float));
	mean_pow_spec = (float *) mxMalloc(((n/2)+1) * sizeof(float));
	twice_inv_noise = (float *) mxMalloc(((n/2)+1) * sizeof(float));
	for (i = 0; i < n; i++) {
		htilde[i] = *(mxGetPr(prhs[5]) + i);
	}
	
	/* call */        
	avg_inv_spec(flo, srate, n, decay, norm, htilde, mean_pow_spec, twice_inv_noise);	

	/* make Matlab output */
	plhs[0] = mxCreateDoubleMatrix((n/2), 1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix((n/2), 1, mxREAL);
	plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
	*mxGetPr(plhs[2]) = *norm;

	for (i = 0; i <= (n/2) ; i++) {
		*(mxGetPr(plhs[0]) + i) = mean_pow_spec[i];
		*(mxGetPr(plhs[1]) + i) = twice_inv_noise[i];
	}
}
 

