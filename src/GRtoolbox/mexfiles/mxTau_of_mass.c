/* copile with:                                                                                                 */
/* mex -I/LIGO/GRASP/include -I/LIGO/frame/v3.60/include mxTau_of_mass.c -L/LIGO/GRASP/lib -lgrasp              */
/*     -L/usr/local/lib -lrecipes_c 			                                                	*/
/*                                                                                                              */
/* example run: [tau0, tau1] = mxTau_of_mass(1.4, 3, pi*100)                                                    */
/*                                                                                                              */
/* Steve Drasco                                                                                                 */
/* Summer 1998                                                                                                  */


#include "mex.h"
#include "grasp.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *m1, *m2, *pf, *tau0, *tau1;	

	/* mex checks */		
	if (nrhs != 3)
		mexErrMsgTxt("Must have three inputs (m1, m2, pf)");
	if (nlhs != 2)
		mexErrMsgTxt("Must have two outputs (tau0, tau1)");
	
	/* arange input and output */
	m1 = mxGetPr(prhs[0]);
	m2 = mxGetPr(prhs[1]);
	pf = mxGetPr(prhs[2]);
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
	tau0 = mxGetPr(plhs[0]);
	tau1 = mxGetPr(plhs[1]);

	/* call GRASP function */	
	tau_of_mass(*m1, *m2, *pf, tau0, tau1);
}
