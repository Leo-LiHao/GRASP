/* compile: mex -g -I/LIGO/GRASP/include  mxDetector_site.c -L/LIGO/GRASP/lib -lgrasp -L/usr/local/lib -lrecipes_c     	*/
/*															*/
/* example run: 		 							                                */
/*															*/
/* Steve Drasco 													*/
/* Summer 1998														*/

#include "mex.h"
#include "grasp.h"
 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
        float   site_parameters[9];
        int     i, site_choice, status, buflength;
	char	*detectors_file;
	char	noise_file[128],whiten_file[128],site_name[128];

	/* mex checks */
        if (nrhs != 2)
                mexErrMsgTxt("Need 2 inputs (detectors_file, site_choice)");
        if (nlhs != 4)
                mexErrMsgTxt("Need 4 outputs (site_parameters, site_name, noise_file, whiten_file)");

	/* build input arguments */
	site_choice = (int) mxGetScalar(prhs[1]);
	buflength = (mxGetN(prhs[0]) * mxGetM(prhs[0])) + 1;
	detectors_file = mxCalloc(buflength, sizeof(char));
	status = mxGetString(prhs[0], detectors_file, buflength);

	/* call detector_site */
	detector_site(detectors_file, site_choice, site_parameters, site_name, noise_file, whiten_file);

	/* build MATLAB output */
	plhs[0] = mxCreateDoubleMatrix(9, 1, mxREAL);
	for (i = 0; i < 9; i++) 
		*(mxGetPr(plhs[0]) + i) = site_parameters[i];
	plhs[1] = mxCreateString(site_name);
	plhs[2] = mxCreateString(noise_file);
	plhs[3] = mxCreateString(whiten_file);

}

