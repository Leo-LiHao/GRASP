/* compile: mex mxSp_filters.c -I/LIGO/GRASP/include -L/LIGO/GRASP/lib -lgrasp -L/usr/local/lib -lrecipes_c		*/
/*															*/
/* example run: compare_chirps.m 							                                */
/*															*/
/* Steve Drasco 													*/
/* Summer 1998														*/

#include "mex.h"
#include "grasp.h"
 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
        float   m1, m2, *ch1, *ch2, fstart, srate, f_c, t_c;
		int     i, n, order;

        /* mex checks */
        if (nrhs != 8)
                mexErrMsgTxt("Need 8 inputs (m1, m2, fstart, n, srate, t_c, f_c, order)");
        if (nlhs != 2)
                mexErrMsgTxt("Need 2 outputs (ch1, ch2)");

        /* assign pointer to input */
        m1 = (float) mxGetScalar(prhs[0]);
        m2 = (float) mxGetScalar(prhs[1]); 
        fstart = (float) mxGetScalar(prhs[2]);
        n = (int) mxGetScalar(prhs[3]);
        srate = (float) mxGetScalar(prhs[4]);
	t_c = (float) mxGetScalar(prhs[5]);
	f_c = (float) mxGetScalar(prhs[6]);
        order = (int) mxGetScalar(prhs[7]);

	/* allocate memory */
        ch1 = (float *)mxMalloc(sizeof(float) * n);
        ch2 = (float *)mxMalloc(sizeof(float) * n);

	/* call sp_filters */

	sp_filters(m1, m2, ch1, ch2, fstart, n, srate, f_c, t_c, order);

	/* set output */	
	plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
        plhs[1] = mxCreateDoubleMatrix(n, 1, mxREAL);
	for(i = 0; i < n; i++){
		*(mxGetPr(plhs[0]) + i) = (double) ch1[i];
		*(mxGetPr(plhs[1]) + i) = (double) ch2[i];
	}
}
 

