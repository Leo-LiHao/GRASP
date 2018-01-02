/* compile: mex -g -I/LIGO/GRASP/include mxCompute_match.c -L/LIGO/GRASP/lib -lgrasp -L/usr/local/lib -lrecipes_c	*/
/*															*/
/* example run: 		 							                                */
/*															*/
/* Steve Drasco 													*/
/* Summer 1998														*/

#include "mex.h"
#include "grasp.h"
 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
        float   outval, m1, m2, *ch0tilde, *ch90tilde, inverse_distance_scale, *twice_inv_noise, flo, s_n0, s_n90, srate;
        int     i, npoint, err_cd_sprs, order;

        /* mex checks */
        if (nrhs != 13)
                mexErrMsgTxt("Need 13 inputs (m1, m2, ch0tilde, ch90tilde, inverse_discance_scale, *twice_inv_noise, flo, s_n0, s_n90, npoint, srate, err_cd_sprs, order)");
        if (nlhs != 1)
                mexErrMsgTxt("Need 1 outputs (outvalue)");

        /* make GRASP input */
	m1 = (float) mxGetScalar(prhs[0]);
	m2 = (float) mxGetScalar(prhs[1]);
	inverse_distance_scale = (float) mxGetScalar(prhs[4]);
	flo = (float) mxGetScalar(prhs[6]);
	s_n0 = (float) mxGetScalar(prhs[7]);
	s_n90 = (float) mxGetScalar(prhs[8]);
	npoint = (int) mxGetScalar(prhs[9]);
	srate = (float) mxGetScalar(prhs[10]);
	err_cd_sprs = (int) mxGetScalar(prhs[11]);
	order = (int) mxGetScalar(prhs[12]);
	ch0tilde = (float *) mxMalloc(npoint * sizeof(float));
	ch90tilde = (float *) mxMalloc(npoint * sizeof(float));
	twice_inv_noise = (float *) mxMalloc(((npoint/2)+1) * sizeof(float));
	for (i = 0; i < npoint; i++) {
		ch0tilde[i] = *(mxGetPr(prhs[2]) + i);
		ch90tilde[i] = *(mxGetPr(prhs[3]) + i);
	}
	for (i = 0; i <= (npoint/2); i++) {
		twice_inv_noise[i] = *(mxGetPr(prhs[5]) + i);
	}

	/* call */ 
	outval = compute_match(m1, m2, ch0tilde, ch90tilde, inverse_distance_scale, twice_inv_noise, flo, s_n0, s_n90,
			npoint, srate, err_cd_sprs, order);

	/* make Matlab output */	
	*mxGetPr(plhs[0]) = outval;
}
 

