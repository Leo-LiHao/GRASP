/* compile: mex mxPhase_frequency.c -I/LIGO/GRASP/include -L/LIGO/GRASP/lib -lgrasp -L/usr/local/lib -lrecipes_c	*/
/*															*/	
/* example run: phase_evoltn.m												*/
/*															*/
/* Steve Drasco 													*/
/* Summer 1998														*/

#include "mex.h"
#include "grasp.h"
 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray
*prhs[])
{
	float	m1, m2, spin1, spin2, *phaseterms, Initial_Freq, Max_Freq_Rqst, Max_Freq_Actual;
	float	Sample_Time, *phase, *frequency, clscnc_time;
	int	n_phaseterms, steps_alloc, steps_filld, err_cd_sprs, i, chirp_ok;

        /* mex checks */
        if (nrhs != 11)
                mexErrMsgTxt("Need 11 inputs (m1, m2, spin1, spin2, n_phaseterms, phaseterms, Initial_Freq, Max_Freq_Rqst, Sample_Time, steps_alloc, err_cd_sprs)");
        if (nlhs != 5)
                mexErrMsgTxt("Need 6 outputs (Max_Freq_Actual, phase, frequency, steps_filld, clscnc_time)");

	/* get scalar inputs */	
        m1 = (float) mxGetScalar(prhs[0]);
        m2 = (float) mxGetScalar(prhs[1]);
	spin1 = (float) mxGetScalar(prhs[2]);
	spin2 = (float) mxGetScalar(prhs[3]);
	n_phaseterms = (int) mxGetScalar(prhs[4]);
	Initial_Freq = (float) mxGetScalar(prhs[6]);
	Max_Freq_Rqst = (float) mxGetScalar(prhs[7]);
	Sample_Time = (float) mxGetScalar(prhs[8]);
        err_cd_sprs = (int) mxGetScalar(prhs[10]);
	phaseterms = (float *) mxMalloc(n_phaseterms*sizeof(float));
		
	/* build phaseterms array */
	for (i = 0; i < n_phaseterms; i++) {
		phaseterms[i] = (float) *(mxGetPr(prhs[5]) + i);
	}

	/* allocate memory */
        if (mxIsEmpty(prhs[9])) {
                phase = frequency = NULL;
        } else {
                steps_alloc = (int) mxGetScalar(prhs[9]);
                phase = (float *)mxMalloc(steps_alloc*sizeof(float));
                frequency = (float *)mxMalloc(steps_alloc*sizeof(float));
        }

	/* call phase_frequency */
	chirp_ok = phase_frequency(m1, m2, spin1, spin2, n_phaseterms, phaseterms,
        	Initial_Freq, Max_Freq_Rqst, &Max_Freq_Actual,
        	Sample_Time, &phase, &frequency,
        	&steps_alloc, &steps_filld, &clscnc_time,
        	err_cd_sprs); 

	/* build output arrays */
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(steps_filld, 1, mxREAL);
	plhs[2] = mxCreateDoubleMatrix(steps_filld, 1, mxREAL);
	plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL);
	plhs[4] = mxCreateDoubleMatrix(1, 1, mxREAL);
	
	/* set scalars */
	*mxGetPr(plhs[0]) = (double) Max_Freq_Actual;
	*mxGetPr(plhs[3]) = (double) steps_filld;
	*mxGetPr(plhs[4]) = (double) clscnc_time;

	/* set vector outputs */
	for (i = 0; i < steps_filld; i++) {
		*(mxGetPr(plhs[1]) + i) = (double) phase[i];
		*(mxGetPr(plhs[2]) + i) = (double) frequency[i];
	}
}
 

