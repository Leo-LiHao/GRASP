/* compile: mex -g -I/LIGO/GRASP/include mxMatch_parab.c -L/LIGO/GRASP/lib -lgrasp -L/usr/local/lib -lrecipes_c		*/
/*															*/
/* example run: 		 							                                */
/*															*/
/* Steve Drasco 													*/
/* Summer 1998														*/

#include "mex.h"
#include "grasp.h"
 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
        float   m1ref, m2ref, matchcont, srate, flo, ftau, semimajor, semiminor, theta, mcoef[3], outval;
        int     i, order, buflength, status;
	char	*noisefile;

        /* mex checks */
        if (nrhs != 8)
                mexErrMsgTxt("Need 8 inputs (m1ref, m2ref, mathcont, order, srate, flo, ftau, noisefile)");
        if (nlhs != 5)
                mexErrMsgTxt("Need 5 outputs (semimajor, semiminor, theta, mcoef, outval)");

        /* make GRASP input */
	m1ref = (float) mxGetScalar(prhs[0]);
	m2ref = (float) mxGetScalar(prhs[1]);
	matchcont = (float) mxGetScalar(prhs[2]);
	order = (int) mxGetScalar(prhs[3]);
	srate = (float) mxGetScalar(prhs[4]);
	flo = (float) mxGetScalar(prhs[5]);
	ftau = (float) mxGetScalar(prhs[6]);
	buflength = (mxGetN(prhs[7]) * mxGetM(prhs[7])) + 1;
	noisefile = mxCalloc(buflength, sizeof(char));
	status = mxGetString(prhs[7], noisefile, buflength);

	/* call math_parab */        
	outval = match_parab(m1ref, m2ref, matchcont, order, srate, flo, ftau, noisefile, &semimajor, &semiminor, &theta, mcoef);

	/* make Matlab output */
	*mxGetPr(plhs[0]) = semimajor;
	*mxGetPr(plhs[1]) = semiminor;
	*mxGetPr(plhs[2]) = theta;
	for (i = 0; i < 3; i++) 
		*(mxGetPr(plhs[3]) + i) = mcoef[i];
	*mxGetPr(plhs[4]) = outval;
}
 

