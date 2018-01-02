/* copile with:                                                                                                 */
/* mex -g -I/LIGO/GRASP/include -I/LIGO/frame/v3.60/include getfri.c -L/LIGO/GRASP/lib -lgrasp 	                */
/* -L/usr/local/lib -lrecipes_c -L/LIGO/frame/v3.60/lib -lFrame                                                 */
/*                                                                                                              */
/* example run: getfri('fri.ascii'); 	                                                                        */
/*                                                                                                              */
/* Steve Drasco                                                                                                 */
/* Summer 1998                                                                                                  */

#include "mex.h"
#include "grasp.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{    
        struct  fgetinput       fgetinput;
        struct  fgetoutput      fgetoutput;
        int     i,buflength, status;
	char	*outfile;
	FILE	*fp;
	
        /* mex checks */
        if (nrhs != 1)
                mexErrMsgTxt("Need 1 input (outfile)");
        if (nlhs != 0)
                mexErrMsgTxt("No outputs");

	/* open file */
	buflength = (mxGetN(prhs[0]) * mxGetM(prhs[0])) + 1;
	outfile = mxCalloc(buflength, sizeof(char));
	status = mxGetString(prhs[0], outfile, buflength);
	fp = fopen(outfile,"w");
	if (fp == NULL)
		mexErrMsgTxt("Error while trying to open file.");

        /* build fget i/o structures */
        fgetinput.nchan = 1;
        fgetinput.npoint = 4096; 
        fgetinput.inlock = 0;
        fgetinput.seek = 1;
        fgetinput.calibrate = 1;
        fgetinput.files = framefiles;
        fgetinput.chnames = (char **)malloc(fgetinput.nchan * sizeof(char *));
        fgetinput.locations = NULL; 
        fgetinput.chnames[0] = "IFO_DMRO";
        fgetoutput.npoint = (int *)malloc(fgetinput.nchan * sizeof(int));
        fgetoutput.ratios = (int *)malloc(fgetinput.nchan * sizeof(int));

        /* call fget_ch */
        fget_ch(&fgetoutput, &fgetinput); 

        /* fill outfile and close it */
        for(i=0;i<fgetoutput.frinum;i++)
		fprintf(fp, "%f\n", fgetoutput.fri[i]);
	fclose(fp);
} 
