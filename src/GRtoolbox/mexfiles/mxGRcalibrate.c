/* copile with:													*/
/* mex -g -I/LIGO/GRASP/include -I/LIGO/frame/v3.60/include mxGRcalibrate.c -L/LIGO/GRASP/lib -lgrasp 		*/
/* -L/usr/local/lib -lrecipes_c -L/LIGO/frame/v3.60/lib -lFrame							*/
/*                                                                                                              */
/* example run: testcalib.m											*/
/*														*/
/* Steve Drasco                                                                                                 */
/* Summer 1998                                                                                                  */

#include "grasp.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	float	*fri, *cplx, srate;
	int	i, frinum, num, method, order;

	/* mex checks */
	if (nrhs != 6)
		mexErrMsgTxt("Need six inputs (fri, frinum, num, srate, method, order)");
	if (nlhs != 1)
		mexErrMsgTxt("Need one output (complex)");

	/* build GRASP input output from Matlab input */
	frinum = (int) mxGetScalar(prhs[1]);
	num = (int) mxGetScalar(prhs[2]);
	srate = (float) mxGetScalar(prhs[3]);
	method = (int) mxGetScalar(prhs[4]);
	order = (int) mxGetScalar(prhs[5]);
	cplx = (float *) mxMalloc((num+2)*sizeof(float));
	fri = (float *)mxMalloc(frinum*sizeof(float));
	for (i = 0; i < frinum; i++) {
		 *(fri + i) = (float) *(mxGetPr(prhs[0]) + i);
	}

	/* call GRASP GRcalibrate */
	GRcalibrate(fri, frinum, num, cplx, srate, method, order);
	
	/* build Matlab output from GRASP output */	
	plhs[0] = mxCreateDoubleMatrix(num+2, 1, mxREAL);
	for(i = 0; i < num + 2; i++) {
		*( mxGetPr(plhs[0]) + i) = cplx[i];
	}
	
}	
