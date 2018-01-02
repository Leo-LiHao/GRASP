/* compile: mex -g -I/LIGO/GRASP/include -I/LIGO/frame/v3.60/include mxTime_inject_chirp.c -L/LIGO/GRASP/lib -lgrasp	*/
/* 		-L/usr/local/lib -lrecipes_c -L/LIGO/frame/v3.60/lib -lFrame                                            */ 
/*															*/	
/* example run:		 							                               		*/
/*															*/
/* Steve Drasco 													*/
/* Summer 1998														*/

#include "mex.h"
#include "grasp.h"
 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
        float   c0, c90, invMpc, *chirp0, *chirp90, *data, *response, *work;
        int     i, offset, n;

        /* mex checks */
        if (nrhs != 9)
                mexErrMsgTxt("Need 9 inputs (c0, c90, offset, invMpc, chirp0, chirp90, response,data, n)");
        if (nlhs != 2)
                mexErrMsgTxt("Need  outputs (data, work)");

        /* make GRASP input */
	c0 = (float) mxGetScalar(prhs[0]);
	c90 = (float) mxGetScalar(prhs[1]);
	offset = (int) mxGetScalar(prhs[2]);
	invMpc = (float) mxGetScalar(prhs[3]);
	n = (int) mxGetScalar(prhs[8]);
	chirp0 = (float *) mxMalloc(n * sizeof(float));
	chirp90 = (float *) mxMalloc(n * sizeof(float));
	response = (float *) mxMalloc( (n+1) * sizeof(float));
	data = (float *) mxMalloc(n * sizeof(float));
	for (i = 0; i < n; i++) {
		chirp0[i] = mxGetPr(prhs[4])[i];
		chirp90[i] = mxGetPr(prhs[5])[i];
		data[i] =  mxGetPr(prhs[7])[i];
		response[i] = mxGetPr(prhs[6])[i];
	}
	response[n] = *(mxGetPr(prhs[6]) + n); /* one extra for you */
	work = (float *) mxMalloc(n * sizeof(float));

	/* call time_inject_chirp */
	time_inject_chirp(c0, c90, offset, invMpc, chirp0, chirp90, data, response, work, n);        

	/* make Matlab output */
	plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(n, 1, mxREAL);
	for (i = 0; i < n; i++) {
		mxGetPr(plhs[0])[i] = data[i];
		mxGetPr(plhs[1])[i] = work[i];
	}
}
