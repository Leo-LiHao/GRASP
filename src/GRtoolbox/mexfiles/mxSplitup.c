/* compile: mex -g -I/LIGO/GRASP/include mxSplitup.c -L/LIGO/GRASP/lib -lgrasp -L/usr/local/lib -lrecipes_c		*/
/*															*/	
/* example run: 		 							                                */
/*															*/
/* Steve Drasco 													*/
/* Summer 1998														*/

#include "mex.h"
#include "grasp.h"
 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray
*prhs[])
{
        float   *working, *template, *r, total;
        int     i, n, p, *indices;

        /* mex checks */
        if (nrhs != 6)
                mexErrMsgTxt("Need 6 inputs (working, template, r, n, total, p)");
        if (nlhs != 1)
                mexErrMsgTxt("Need 1 outputs (indices)");

        /* make GRASP input */
	n = (int) mxGetScalar(prhs[3]);
	total = (float) mxGetScalar(prhs[4]);
	p = (int) mxGetScalar(prhs[5]);
	working = (float *) mxMalloc(n * sizeof(float));
	template = (float *) mxMalloc(n * sizeof(float));
	r = (float *) mxMalloc(((n/2)+1) * sizeof(float));
	for (i = 0; i < n; i++) {
		working[i] = *(mxGetPr(prhs[0]) + i);
		template[i] = *(mxGetPr(prhs[1]) + i);
	}
	for (i = 0; i <= (n/2); i++) {
		r[i] = *(mxGetPr(prhs[2]) + i);
	}
	indices = (int *) mxMalloc(p * sizeof(int));

	/* call */        
	splitup(working, template, r, n, total, p, indices);

	/* make Matlab output */
	plhs[0] = mxCreateDoubleMatrix(p, 1, mxREAL);
	for (i = 0; i < p; i++) {
		*(mxGetPr(plhs[0]) +i) = indices[i];
	}
}
 

