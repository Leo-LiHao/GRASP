/* compile: mex -g -I/LIGO/GRASP/include mxTemplate_area.c -L/LIGO/GRASP/lib -lgrasp -L/usr/local/lib -lrecipes_c	*/
/*															*/
/* example run: 		 							                                */
/*															*/
/* Steve Drasco 													*/
/* Summer 1998														*/

#include "mex.h"
#include "grasp.h"
 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	struct	Scope Grid;
	mxArray	*n_tmplt_array, *m_mn_array, *m_mx_array, *theta_array, *dp_array, *dq_array, *f_start_array;
	float	area;
	const	char *grid_fields[8];

        /* mex checks */
        if (nrhs != 1)
                mexErrMsgTxt("Need 1 inputs (Grid - scope structure)");
        if (nlhs != 2)
                mexErrMsgTxt("Need  outputs (Grid - scope structure, area)");

        /* make GRASP input */
	if( !mxIsEmpty(mxGetField(prhs[0], 0, "n_tmplt")) ) 
		Grid.n_tmplt = (int) mxGetScalar( mxGetField(prhs[0], 0, "n_tmplt") );
        if( !mxIsEmpty(mxGetField(prhs[0], 0, "m_mn")) ) 
		Grid.m_mn = (float) mxGetScalar( mxGetField(prhs[0], 0, "m_mn") );
        if( !mxIsEmpty(mxGetField(prhs[0], 0, "m_mx") ) )
		Grid.m_mx = (float) mxGetScalar( mxGetField(prhs[0], 0, "m_mx") );
        if( !mxIsEmpty(mxGetField(prhs[0], 0, "theta")) )
		Grid.theta = (float) mxGetScalar( mxGetField(prhs[0], 0, "theta") );
        if( !mxIsEmpty(mxGetField(prhs[0], 0, "dp")) )
		Grid.dp = (float) mxGetScalar( mxGetField(prhs[0], 0, "dp") );
        if( !mxIsEmpty(mxGetField(prhs[0], 0, "dq")) )
		Grid.dq = (float) mxGetScalar( mxGetField(prhs[0], 0, "dq") );
        if( !mxIsEmpty(mxGetField(prhs[0], 0, "f_start")) )
		Grid.f_start = (float) mxGetScalar( mxGetField(prhs[0], 0, "f_start") );
	
	/* call template_area */        
	area = template_area(&Grid);

	/* set output value */
	plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
	*(mxGetPr(plhs[1])) = area;

	/* field names */	
	grid_fields[0] = "n_tmplt";
        grid_fields[1] = "m_mn";
        grid_fields[2] = "m_mx";
        grid_fields[3] = "theta";
        grid_fields[4] = "dp";
        grid_fields[5] = "dq";
        grid_fields[6] = "f_start";
        grid_fields[7] = "templates";

	plhs[0] = mxCreateStructMatrix(1, 1, 8, grid_fields);
	
	n_tmplt_array = mxCreateDoubleMatrix(1,1,mxREAL);
	m_mn_array = mxCreateDoubleMatrix(1,1,mxREAL);
	m_mx_array = mxCreateDoubleMatrix(1,1,mxREAL);
	theta_array = mxCreateDoubleMatrix(1,1,mxREAL);
	dp_array = mxCreateDoubleMatrix(1,1,mxREAL);
	dq_array = mxCreateDoubleMatrix(1,1,mxREAL);
	f_start_array = mxCreateDoubleMatrix(1,1,mxREAL);
        
	*mxGetPr(n_tmplt_array) =  Grid.n_tmplt;
	*mxGetPr(m_mn_array) =  Grid.m_mn;
	*mxGetPr(m_mx_array) =  Grid.m_mx;
	*mxGetPr(theta_array) =  Grid.theta;
	*mxGetPr(dp_array) =  Grid.dp;
	*mxGetPr(dq_array) =  Grid.dq;
	*mxGetPr(f_start_array) =  Grid.f_start;

	mxSetField(plhs[0], 0, "n_tmplt", n_tmplt_array);
	mxSetField(plhs[0], 0, "m_mn", m_mn_array);
	mxSetField(plhs[0], 0, "m_mx", m_mx_array);
	mxSetField(plhs[0], 0, "theta", theta_array);
	mxSetField(plhs[0], 0, "dp", dp_array);
	mxSetField(plhs[0], 0, "dq", dq_array);
	mxSetField(plhs[0], 0, "f_start", f_start_array);
	mxSetField( plhs[0], 0, "templates", mxDuplicateArray(mxGetField(prhs[0], 0, "templates")) );
}
