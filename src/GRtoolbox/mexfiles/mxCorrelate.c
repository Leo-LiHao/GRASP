/* compile: mex -g -I/LIGO/GRASP/include mxCorrelate.c -L/LIGO/GRASP/lib -lgrasp -L/usr/local/lib -lrecipes_c		*/
/*															*/	
/* example run: 		 							                                */
/*															*/
/* Steve Drasco 													*/
/* Summer 1998														*/

#include "mex.h"
#include "grasp.h"
 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
        float	*s, *h, *c, *r;
        int	i, n;

        /* mex checks */
        if (nrhs != 4)
                mexErrMsgTxt("Need 4 inputs (h, c, r, n)");
        if (nlhs != 1)
                mexErrMsgTxt("Need 1 outputs (s)");

        /* make GRASP input */
	n = (int) mxGetScalar(prhs[3]);
	c = (float *) mxMalloc(n * sizeof(float));
	h = (float *) mxMalloc(n * sizeof(float));
	r = (float *) mxMalloc(((n/2)+1) * sizeof(float));
	s = (float *) mxMalloc(n * sizeof(float));
	for (i = 0; i < n; i++) {
		h[i] = *(mxGetPr(prhs[0]) + i);
		c[i] = *(mxGetPr(prhs[1]) + i);
	}
	for (i = 0; i <= n/2; i++) {
		r[i] = *(mxGetPr(prhs[2]) + i);
	}

	/* call */        
	correlate(s, h, c, r, n);

	/* make Matlab output */
	plhs[0] = mxCreateDoubleMatrix(2, 1, mxREAL);
	for (i = 0; i < n; i++) {
                *(mxGetPr(plhs[0]) + i) = s[i];
        }	
}
 

