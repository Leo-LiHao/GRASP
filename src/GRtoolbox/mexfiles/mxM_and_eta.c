/* copile with:                                                                                                 */
/* mex -I/LIGO/GRASP/include mxM_and_eta.c -L/LIGO/GRASP/lib -lgrasp -L/usr/local/lib -lrecipes_c		*/
/*                                                                                                              */
/* example run: [tau0, tau1] = mxTau_of_mass(1.4, 3, pi*100);							*/
/*		[m, eta] = mxM_and_eta(tau0, tau1, 0.5, 3.0, pi*60)						*/
/*	[ NOTE : THESE RESULTS GIVE M = ETA = 0 ... TRY TO FIND SOME THAT DON'T ]				*/
/*                                                                                                              */
/* Steve Drasco                                                                                                 */
/* Summer 1998                                                                                                  */


#include "mex.h"
#include "grasp.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *Mmin, *Mmax, *pf, *tau0, *tau1, *eta, *M;	

	/* mex checks */		
	if (nrhs != 5)
		mexErrMsgTxt("Must have 5 inputs (tau0, tau1, Mmin, Mmax, pf)");
	if (nlhs != 2)
		mexErrMsgTxt("Must have two outputs (m, eta)");
	
	/* arange input and output */
	tau0 = mxGetPr(prhs[0]);
	tau1 = mxGetPr(prhs[1]);
	Mmin = mxGetPr(prhs[2]);
	Mmax = mxGetPr(prhs[3]);
	pf = mxGetPr(prhs[4]);
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
	M = mxGetPr(plhs[0]);
	eta = mxGetPr(plhs[1]);
	
	/* call GRASP m_and_eta */	
	m_and_eta(*tau0, *tau1, M, eta, *Mmin, *Mmax, *pf);
}
