/* compile:  mex mxMake_filters.c -I/LIGO/GRASP/include -L/LIGO/GRASP/lib -lgrasp -L/usr/local/lib -lrecipes_c		*/
/*															*/
/* example run: [ch1, ch2, filled, t_coal]=mxMake_filters(10,10,40,32768,10000,0,4)					*/
/*															*/
/* Steve Drasco 													*/
/* Summer 1998														*/

#include "mex.h"
#include "grasp.h"
 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
        float   m1, m2, fstart, srate, t_coal, *ch1, *ch2;
        int     filled, i, n, err_cd_sprs, order;

        /* mex checks */
        if (nrhs != 7)
                mexErrMsgTxt("Need 7 inputs (m1, m2, fstart, length, srate, err_cd_sprs, order)");
        if (nlhs != 4)
                mexErrMsgTxt("Need 4 outputs (ch1, ch2, filled, t_coal)");

        /* assign pointer to input */
	m1 = (float) mxGetScalar(prhs[0]);
        m2 = (float) mxGetScalar(prhs[1]);
	fstart = (float) mxGetScalar(prhs[2]);
	n = (int) mxGetScalar(prhs[3]);
	srate = (float) mxGetScalar(prhs[4]);
        err_cd_sprs = (int) mxGetScalar(prhs[5]);
        order = (int) mxGetScalar(prhs[6]);

	/* allocate memory */
        ch1 = (float *)mxMalloc(sizeof(float) * n);
        ch2 = (float *)mxMalloc(sizeof(float) * n);
        
	make_filters(m1, m2, ch1, ch2, fstart, n, srate, &filled, &t_coal, err_cd_sprs, order);

	/* set scalar output */ 
        plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
	plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL);
        *mxGetPr(plhs[2]) = (double) filled;
	*mxGetPr(plhs[3]) = (double) t_coal;

        plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
        plhs[1] = mxCreateDoubleMatrix(n, 1, mxREAL);
        for(i = 0; i < n; i++){
                *(mxGetPr(plhs[0]) + i) = (double) ch1[i]; 
                *(mxGetPr(plhs[1]) + i) = (double) ch2[i]; 
        }
}
 

