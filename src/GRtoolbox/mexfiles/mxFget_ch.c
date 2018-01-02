/* copile with:													*/
/* mex -g -I/LIGO/GRASP/include -I/LIGO/frame/v3.60/include mxFget_ch.c -L/LIGO/GRASP/lib -lgrasp 		*/
/* -L/usr/local/lib -lrecipes_c -L/LIGO/frame/v3.60/lib -lFrame							*/
/*                                                                                                              */
/* example run: testget												*/
/*														*/
/* Steve Drasco                                                                                                 */
/* Summer 1998                                                                                                  */

#include "grasp.h"
#include "mex.h"
#include "string.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	mxArray	*chnames_field_ptr, *cell_array_ptr;
	mxArray	*tstart_array, *srate_array, *npoint_array, *ratios_array, *discarded_array;
	mxArray *tfirst_array, *dt_array, *lostlock_array, *lastlock_array, *returnval_array;
	mxArray *frinum_array, *fri_array, *tcalibrate_array, *locklow_array, *lockhi_array;
	mxArray *data_array, *cell_ptr, *data_matrix;
	struct	fgetinput	fgetinput;
	struct	fgetoutput	fgetoutput;
	int	i, j, dims[2], buflength, status;
	const	char	*foutnames[16];

	/* mex checks */
	if (nrhs != 1)
		mexErrMsgTxt("Need one  input");
	else if(!mxIsStruct(prhs[0]))
		mexErrMsgTxt("Input must be a structure");
	if (nlhs != 1)
		mexErrMsgTxt("Need one output");
	
	/* build GRASP fgetinput from Matlab fgetinput */
	fgetinput.nchan = (int) *mxGetPr( mxGetField(prhs[0], 0, "nchan") );
	fgetinput.npoint = (int) *mxGetPr( mxGetField(prhs[0], 0, "npoint") );
	fgetinput.inlock = (int) *mxGetPr( mxGetField(prhs[0], 0, "inlock") );
	fgetinput.seek = (int) *mxGetPr( mxGetField(prhs[0], 0, "seek" ) );
	fgetinput.calibrate = (int) *mxGetPr( mxGetField(prhs[0], 0, "calibrate") );
	fgetinput.files = framefiles;
	fgetinput.chnames = (char **)mxMalloc(fgetinput.nchan * sizeof(char *));
	if (fgetinput.seek == 0) {
		fgetinput.locations = (short **)mxMalloc(fgetinput.nchan * sizeof(short *));	
		for(i = 0; i < fgetinput.nchan; i++)
			fgetinput.locations[i] = (short *)mxMalloc(fgetinput.npoint * sizeof(short));
	} else {
		fgetinput.locations = NULL;
	}
	chnames_field_ptr = mxGetField(prhs[0],0, "chnames" );
	for(i = 0; i < (mxGetN(chnames_field_ptr) * mxGetM(chnames_field_ptr)); i++) {
		cell_array_ptr = mxGetCell(chnames_field_ptr, i);
		buflength = (mxGetN(cell_array_ptr) * mxGetM(cell_array_ptr)) + 1;
		fgetinput.chnames[i] = mxCalloc(buflength, sizeof(char));
		status = mxGetString(cell_array_ptr, fgetinput.chnames[i], buflength);
	}

	/* build GRASP fgetoutput from Matlab fgetoutput */
	fgetoutput.npoint = (int *)mxMalloc(fgetinput.nchan * sizeof(int));
	fgetoutput.ratios = (int *)mxMalloc(fgetinput.nchan * sizeof(int));
		
	/* call GRASP fget_ch */
	fget_ch(&fgetoutput, &fgetinput); 
	
	/* TO CHECK THE REAL VALUES *********************************************	
	mexPrintf("\n|------ These values come from GRASP itself--------|\n");
	mexPrintf("\ntstart = %f\nsrate = %f\ndiscarded = %d\ntfirst = %f\ndt = %f\nlostlock = %f\n",
                fgetoutput.tstart, fgetoutput.srate, fgetoutput.discarded,
                fgetoutput.tfirst, fgetoutput.dt, fgetoutput.lostlock);
	mexPrintf("lastlock = %f \nreturnval = %d\nfrinum = %d\ntcalibrate  = %d\nlocklow = %d\nlockhi = %d\n",
		fgetoutput.lastlock, fgetoutput.returnval, fgetoutput.frinum,
                fgetoutput.tcalibrate, fgetoutput.locklow, fgetoutput.lockhi);
	mexPrintf("\n|------ Next values come mex file --------|\n");
	************************************************************************* */
	

	/* build Matlab fgetoutput from GRASP fgetoutput */
	foutnames[0]="tstart";
	foutnames[1]="srate";
	foutnames[2]="npoint";
	foutnames[3]="ratios";
	foutnames[4]="discarded";
	foutnames[5]="tfirst";
	foutnames[6]="dt";
	foutnames[7]="lostlock";
	foutnames[8]="lastlock";
	foutnames[9]="returnval";
	foutnames[10]="frinum";
	foutnames[11]="fri";
	foutnames[12]="tcalibrate";
	foutnames[13]="locklow";
	foutnames[14]="lockhi";
	foutnames[15]="data";
	plhs[0] = mxCreateStructMatrix(1, 1, 16, foutnames);
	
	tstart_array = mxCreateDoubleMatrix(1,1,mxREAL);
	srate_array = mxCreateDoubleMatrix(1,1,mxREAL);
	npoint_array = mxCreateDoubleMatrix(1,1,mxREAL);
	ratios_array = mxCreateDoubleMatrix(1,1,mxREAL);
	discarded_array = mxCreateDoubleMatrix(1,1,mxREAL);
	tfirst_array = mxCreateDoubleMatrix(1,1,mxREAL);
	dt_array = mxCreateDoubleMatrix(1,1,mxREAL);
	lostlock_array = mxCreateDoubleMatrix(1,1,mxREAL);
	lastlock_array = mxCreateDoubleMatrix(1,1,mxREAL);
	returnval_array = mxCreateDoubleMatrix(1,1,mxREAL);
	frinum_array = mxCreateDoubleMatrix(1,1,mxREAL);
	fri_array = mxCreateDoubleMatrix(1,1,mxREAL);
	tcalibrate_array = mxCreateDoubleMatrix(1,1,mxREAL);
	locklow_array = mxCreateDoubleMatrix(1,1,mxREAL);
	lockhi_array = mxCreateDoubleMatrix(1,1,mxREAL);

	*mxGetPr(tstart_array) =  fgetoutput.tstart;
	*mxGetPr(srate_array) = fgetoutput.srate;
	*mxGetPr(discarded_array) = fgetoutput.discarded;
	*mxGetPr(tfirst_array) = fgetoutput.tfirst;
	*mxGetPr(dt_array) = fgetoutput.dt;
	*mxGetPr(lostlock_array) = fgetoutput.lostlock;
	*mxGetPr(lastlock_array) = fgetoutput.lastlock;
	*mxGetPr(returnval_array) = fgetoutput.returnval;
	*mxGetPr(frinum_array) = fgetoutput.frinum;
	*mxGetPr(tcalibrate_array) = fgetoutput.tcalibrate;
	*mxGetPr(locklow_array) = fgetoutput.locklow;
	*mxGetPr(lockhi_array) = fgetoutput.lockhi;			 

	/* scalar fields */
	mxSetField(plhs[0], 0, "tstart", tstart_array);	
	mxSetField(plhs[0], 0, "srate",srate_array);
	mxSetField(plhs[0], 0, "discarded",discarded_array);
	mxSetField(plhs[0], 0, "tfirst",tfirst_array);
	mxSetField(plhs[0], 0, "dt",dt_array);
	mxSetField(plhs[0], 0, "lostlock",lostlock_array);
	mxSetField(plhs[0], 0, "lastlock",lastlock_array);
	mxSetField(plhs[0], 0, "returnval",returnval_array);
	mxSetField(plhs[0], 0, "frinum",frinum_array);
	mxSetField(plhs[0], 0, "tcalibrate",tcalibrate_array);
	mxSetField(plhs[0], 0, "locklow",locklow_array);
	mxSetField(plhs[0], 0, "lockhi",lockhi_array);

	/* data field */		
        if (fgetinput.seek == 0) {
                dims[0] = 1;
                dims[1] = fgetinput.nchan;
                data_array = mxCreateCellArray(2, dims);
                for(i = 0; i < fgetinput.nchan; i++) {
                        cell_ptr = mxGetCell(data_array, i);
                        data_matrix = mxCreateDoubleMatrix(fgetoutput.npoint[i], 1, mxREAL);

                        for(j = 0; j < fgetoutput.npoint[i]; j++)  {
                                 *( mxGetPr(data_matrix) + j) = *(fgetinput.locations[i] + j);
                        }

                        mxSetCell(data_array, i, data_matrix);
                }
		mxSetField(plhs[0], 0, "data", data_array);
        }
	
	/* ratios field */
        ratios_array = mxCreateDoubleMatrix(fgetinput.nchan, 1, mxREAL);
        for(i = 0; i < fgetinput.nchan; i++)
              	*( mxGetPr(ratios_array) + i) = fgetoutput.ratios[i];
        mxSetField(plhs[0], 0, "ratios", ratios_array);

	/* npoint field */
	npoint_array = mxCreateDoubleMatrix(fgetinput.nchan, 1, mxREAL);
	for(i = 0; i < fgetinput.nchan; i++)
                *( mxGetPr(npoint_array) + i) = fgetoutput.npoint[i];
        mxSetField(plhs[0], 0, "npoint", npoint_array);
	
	/* fri field */
	fri_array = mxCreateDoubleMatrix(fgetoutput.frinum, 1, mxREAL);
        for(i = 0; i < fgetoutput.frinum; i++)
                *( mxGetPr(fri_array) + i) = fgetoutput.fri[i];
        mxSetField(plhs[0], 0, "fri", fri_array);
}	
