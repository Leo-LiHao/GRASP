/* copile with:													*/
/* mex -g -I/LIGO/GRASP/include -I/LIGO/frame/v3.60/include mxGRnormalize.c -L/LIGO/GRASP/lib -lgrasp 		*/
/* -L/usr/local/lib -lrecipes_c -L/LIGO/frame/v3.60/lib -lFrame							*/
/*                                                                                                              */
/* example run: power_spectrumF.m										*/
/*														*/
/* Steve Drasco                                                                                                 */
/* Summer 1998                                                                                                  */

#include "grasp.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	float	*fri, *response, srate, npoint;
	int	i, frinum, num, method, order;

	/* mexchecks */
	if (nrhs != 4)
		mexErrMsgTxt("Need four inputs (fri, frinum, npoint, srate)");
	if (nlhs != 1)
		mexErrMsgTxt("Need one output (response)");

	/* build GRASP input */
	frinum = (int) mxGetScalar(prhs[1]);
	fri = (float *)mxMalloc(frinum*sizeof(float));
	for (i = 0; i < frinum; i++) {
                 *(fri + i) = (float) *(mxGetPr(prhs[0]) + i);
        }
	npoint = (int) mxGetScalar(prhs[2]);
	srate = (float) mxGetScalar(prhs[3]);
	response = (float *) mxMalloc((npoint + 2)*sizeof(float));

	/* call GRASP GRnormalize */
	GRnormalize(fri, frinum, npoint, srate, response);
	

	/* build Matlab output */
	plhs[0] = mxCreateDoubleMatrix(npoint+2, 1, mxREAL);
	for(i = 0; i < npoint + 2; i++) {
		mxGetPr(plhs[0])[i] = response[i];
	} 
}	
