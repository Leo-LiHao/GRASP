/* compile: mex -g -I/LIGO/GRASP/include mxOrthonormalize.c -L/LIGO/GRASP/lib -lgrasp -L/usr/local/lib -lrecipes_c	*/
/*															*/	
/* example run: 		 							                                */
/*															*/
/* Steve Drasco 													*/
/* Summer 1998														*/

#include "mex.h"
#include "grasp.h"
 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
        float   *ch0tilde, *ch90tilde, *twice_inv_noise, *n0, *n90;
        int     i, n;

        /* mex checks */
        if (nrhs != 4)
                mexErrMsgTxt("Need 4 inputs (ch0tilde, ch90tilde, twice_inv_noise, n)");
        if (nlhs != 3)
                mexErrMsgTxt("Need 3 outputs (n0, n90, ch90tilde)");

        /* make GRASP input */
	n = (int) mxGetScalar(prhs[3]);
	ch0tilde = (float *) mxMalloc(n * sizeof(float));
	ch90tilde = (float *) mxMalloc(n * sizeof(float));
	n0 = (float *) mxMalloc(n * sizeof(float));
	n90 = (float *) mxMalloc(n * sizeof(float));
	twice_inv_noise = (float *) mxMalloc(((n/2)+1) * sizeof(float));
	for ( i = 0; i < n; i++) {
		ch0tilde[i] = *(mxGetPr(prhs[0]) + i);
		ch90tilde[i] = *(mxGetPr(prhs[1]) + i);
	}
	for (i = 0; i <= (n/2); i++) {
		twice_inv_noise[i] = *(mxGetPr(prhs[2]) + i);
	}
		

	/* call */        
	orthonormalize(ch0tilde, ch90tilde, twice_inv_noise, n, n0, n90);

	/* make Matlab output */
	plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(n, 1, mxREAL);
	plhs[2] = mxCreateDoubleMatrix((n/2)+1, 1, mxREAL);
	for (i = 0; i < n; i++) {
               *(mxGetPr(plhs[0]) + i) = n0[i];
               *(mxGetPr(plhs[1]) + i) = n90[i];
        }
	for (i = 0; i <= (n/2); i++) {
		*(mxGetPr(plhs[2]) + i) = ch90tilde[i];
        }
}
