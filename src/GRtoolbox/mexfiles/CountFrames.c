/* CountFrames.c  
 *	This programs counts the number of frames in a frame file. 
 *	It prints the Frame numbers of the first and last frames.
 *
 *      input: file path as string
 *	
 *	output: number of frames
 *
 * Steve Drasco, Summer 1999
 */

#include "FrameL.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int     buflength, status;
	double	nFrame;
	char	*FrameFile;

        /* mex checks */
        if (nrhs != 1) mexErrMsgTxt("Need 1 input (location of frame file)");
        if (nlhs != 1) mexErrMsgTxt("Need 1 output (number of frames in file)");

	/* read in the file path */
	buflength = (mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;
        FrameFile = (char *) mxCalloc(buflength, sizeof(char));
        if( (status = mxGetString(prhs[0], FrameFile, buflength)) != 0) mexErrMsgTxt("Could not convert string data.");

	/* call the counting routine */
	if ( (nFrame = (double) CountFrames(FrameFile)) == 0.) mexErrMsgTxt("");

	/* arrange output */
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	*( mxGetPr(plhs[0]) ) = nFrame; 
}

int CountFrames(char *FrameFile)
{
	int    nFrame=0;
	struct FrFile *iFile;
	struct FrameH *frame;

	/* open the frame file */	
	if ( (iFile = FrFileINew(FrameFile)) == NULL ) {
		mexPrintf("Can not open file: %s\n",FrameFile);
		return 0;
	}

	/* Count the frames */
	while( (frame=FrameRead(iFile)) != NULL ) nFrame++;

	/* Close the file and exit */
	FrFileIEnd(iFile);
	return nFrame;
}
