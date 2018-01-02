/* compile: mex -g -I/LIGO/GRASP/include -I/LIGO/frame/v3.60/include mxFreq_inject_chirp.c -L/LIGO/GRASP/lib -lgrasp	*/
/* 		-L/usr/local/lib -lrecipes_c -L/LIGO/frame/v3.60/lib -lFrame                                            */ 
/*															*/	
/* example run: 		 							                                */
/*															*/
/* Steve Drasco 													*/
/* Summer 1998														*/

#include "mex.h"
#include "grasp.h"
 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
        float   c0, c90, invMpc, *ch0tilde, *ch90tilde, *htilde;
        int     i, offset, n;

        /* mex checks */
        if (nrhs != 7)
                mexErrMsgTxt("Need 7 inputs (c0, c90, offset, invMpc, ch0tilde, ch90tilde, n)");
        if (nlhs != 1)
                mexErrMsgTxt("Need 1 outputs (htilde)");

        /* make GRASP input */
	c0 = (float) mxGetScalar(prhs[0]);
	c90 = (float) mxGetScalar(prhs[1]);
	offset = (int) mxGetScalar(prhs[2]);
	invMpc = (float) mxGetScalar(prhs[3]);
	n = (int) mxGetScalar(prhs[6]);
	ch0tilde = (float *) mxMalloc(n * sizeof(float));
	ch90tilde = (float *) mxMalloc(n * sizeof(float));
	for (i = 0; i < n; i++) {
		ch0tilde[i] = *(mxGetPr(prhs[4]) +i);
		ch90tilde[i] = *(mxGetPr(prhs[5]) +i);
	}
	htilde = (float *) mxMalloc(n * sizeof(float));

	/* call */        
	freq_inject_chirp(c0, c90, offset, invMpc, ch0tilde, ch90tilde, htilde, n);

	/* make Matlab output */	
	plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
	for (i = 0; i < n; i++) {
		*(mxGetPr(plhs[0]) + i) = htilde[i];
	}
}
 

