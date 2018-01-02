/* compile: mex -g -I/LIGO/GRASP/include mxSplitup_freq2.c -L/LIGO/GRASP/lib -lgrasp -L/usr/local/lib -lrecipes_c	*/
/*															*/	
/* example run: 		 							                                */
/*															*/
/* Steve Drasco 													*/
/* Summer 1998														*/

#include "mex.h"
#include "grasp.h"
 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	float	c0, c90, *chirp0, *chirp90, norm, *twice_inv_noise, *stats, *working, *htilde, rsquared;
        int     i, n, offset, p, *indices;

        /* mex checks */
        if (nrhs != 9)
                mexErrMsgTxt("Need 9 inputs (c0, c90, chirp0, chirp90, norm, twice_inv_noise, n, offset, p)");
        if (nlhs != 5)
                mexErrMsgTxt("Need 5 outputs (indices, stats, working, htilde, rsquared)");

        /* make GRASP input */
	c0 = (float) mxGetScalar(prhs[0]);
	c90 = (float) mxGetScalar(prhs[1]);
	norm = (float) mxGetScalar(prhs[4]);
	n = (int) mxGetScalar(prhs[6]);
	offset = (int) mxGetScalar(prhs[7]);
	p = (int) mxGetScalar(prhs[8]);
	chirp0 = (float *) mxMalloc(n * sizeof(float));
	chirp90 = (float *) mxMalloc(n * sizeof(float));
	twice_inv_noise = (float *) mxMalloc(((n/2)+1) * sizeof(float));
	for (i = 0; i < n; i++) {
		chirp0[i] = *(mxGetPr(prhs[2]) + i);
		chirp90[i] = *(mxGetPr(prhs[3]) + i);
	}
	for (i = 0; i <= (n/2); i++) {
		twice_inv_noise[i] = *(mxGetPr(prhs[5]) + i);
	}
	indices = (int *) mxMalloc(p * sizeof(int));
	stats = (float *) mxMalloc(p * sizeof(float));
	working = (float *) mxMalloc(n * sizeof(float));
	htilde = (float *) mxMalloc(n * sizeof(float));

	/* call splitup_freq2 */        
	rsquared = splitup_freq2(c0, c90, chirp0, chirp90, norm, twice_inv_noise, n, offset, p, indices, stats, working, htilde);

	/* make Matlab output */
	plhs[0] = mxCreateDoubleMatrix(p, 1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(p, 1, mxREAL);
	plhs[2] = mxCreateDoubleMatrix(n, 1, mxREAL);
	plhs[3] = mxCreateDoubleMatrix(n, 1, mxREAL);
	plhs[4] = mxCreateDoubleMatrix(1, 1, mxREAL);
	*mxGetPr(plhs[4]) = rsquared;

	for (i = 0; i < p; i++) {
		*(mxGetPr(plhs[0]) +i) = indices[i];
		*(mxGetPr(plhs[1]) +i) = stats[i];
	}
	for (i = 0; i < n; i++) {
                *(mxGetPr(plhs[2]) +i) = working[i];
                *(mxGetPr(plhs[3]) +i) = htilde[i];
        }
}
 

