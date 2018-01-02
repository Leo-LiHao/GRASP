/* compile: mex -g -I/LIGO/GRASP/include -I/LIGO/frame/v3.60/include mxFind_chirp.c -L/LIGO/GRASP/lib -lgrasp		*/
/* 		-L/usr/local/lib -lrecipes_c -L/LIGO/frame/v3.60/lib -lFrame                                            */ 
/*															*/
/* example run: .m		 							                                */
/*															*/
/* Steve Drasco 													*/
/* Summer 1998														*/

#include "mex.h"
#include "grasp.h"
 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
        float   *htilde, *ch0tilde, *ch90tilde, *twice_inv_noise, n0, n90, *output0, *output90, snr_max, c0, c90, var;
        int     i, n, chirplen, offset;

        /* mex checks */
        if (nrhs != 8)
                mexErrMsgTxt("Need 8 inputs (htilde, ch0tilde, ch90tilde, twice_inv_noise, n0, n90, n, chirplen)");
        if (nlhs != 7)
                mexErrMsgTxt("Need 7 outputs (output0, output90, offset, snr_max, c0, c90, var)");

        /* make GRASP input */
	n0 = (float) mxGetScalar(prhs[4]);
	n90 = (float) mxGetScalar(prhs[5]);
	n = (int) mxGetScalar(prhs[6]);
	chirplen = (float) mxGetScalar(prhs[7]);
	htilde = (float *) mxMalloc(n * sizeof(float));
	ch0tilde = (float *) mxMalloc(n * sizeof(float));
	ch90tilde = (float *) mxMalloc(n * sizeof(float));
	output0 = (float *) mxMalloc(n * sizeof(float));
	output90 = (float *) mxMalloc(n * sizeof(float));
	twice_inv_noise = (float *) mxMalloc( ((n/2)+1) * sizeof(float));
	for (i = 0; i < n; i++) {
		ch0tilde[i] = *(mxGetPr(prhs[1]) + i);
		ch90tilde[i] = *(mxGetPr(prhs[2]) + i);
		htilde[i] = *(mxGetPr(prhs[0]) + i);
	}
	for (i = 0; i <= (n/2); i++) {
		twice_inv_noise[i] = *(mxGetPr(prhs[3]) + i);
	}

	/* call find_chirp */        
	findchirp(htilde, ch0tilde, ch90tilde, twice_inv_noise, n0, n90, output0, output90,
		 n, chirplen, &offset, &snr_max, &c0, &c90, &var);
	
	/* make Matlab output */
	plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(n, 1, mxREAL);
	plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
	plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL);
	plhs[4] = mxCreateDoubleMatrix(1, 1, mxREAL);
	plhs[5] = mxCreateDoubleMatrix(1, 1, mxREAL);
	plhs[6] = mxCreateDoubleMatrix(1, 1, mxREAL);

	*mxGetPr(plhs[2]) = offset ;
	*mxGetPr(plhs[3]) = snr_max;
	*mxGetPr(plhs[4]) = c0;
	*mxGetPr(plhs[5]) = c90;
	*mxGetPr(plhs[6]) = var;

	for (i = 0; i < n; i++) {
		*(mxGetPr(plhs[0]) + i) = output0[i];
		*(mxGetPr(plhs[1]) + i) = output90[i];
	}

}
 

