/* compile: mex -g -I/LIGO/GRASP/include mxGet_cubic.c -L/LIGO/GRASP/lib -lgrasp -L/usr/local/lib -lrecipes_c		*/
/*															*/
/* example run: 		 							                                */
/*															*/
/* Steve Drasco 													*/
/* Summer 1998														*/

#include "mex.h"
#include "grasp.h"
 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	struct	cubic_grid grid;
	float	m1, m2, coef[10];
        int	i, j, k, dims[3], returnval;

	/* mex checks */
        if (nrhs != 3)
                mexErrMsgTxt("Need 3 input (m1, m2, grid)");
        if (nlhs != 2)
                mexErrMsgTxt("Need 2 outputs (coef, returnval)");

	/* get scalars */
	m1 = mxGetScalar(prhs[0]);
	m2 = mxGetScalar(prhs[1]);
		
        /* allocate memory for grid.coef */
        dims[0] = dims[1] = grid.n;
        dims[2] = 10;
        grid.coef = (float ***) mxMalloc(dims[2]*sizeof(float **));
        for (i = 0; i < dims[2]; i++) {
                grid.coef[i] = (float **) mxMalloc(dims[1]*sizeof(float *));
                for (j = 0; j < dims[1]; j++) {
                        grid.coef[i][j] = (float *) mxMalloc(dims[0]*sizeof(float));
                }
        }
        /* fill memory of grid.coef */
        for (i = 0; i < dims[0]; i++) {
                for (j = 0; j < dims[1]; j++) {
                        for (k = 0; k < dims[2]; k++) {
                                grid.coef[k][j][i] = *(mxGetPr(prhs[0]) + (dims[2]*dims[1]*i + dims[2]*j + k) );
                        }
                }
        }
	
	/* call get_cubic */
	returnval = get_cubic(m1, m2, grid, coef);

	/* make output */	
	plhs[0] = mxCreateDoubleMatrix(10, 1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
	*mxGetPr(plhs[1]) = returnval;
        for (i = 0; i < dims[0]; i++) {
		*(mxGetPr(prhs[0]) + i) = coef[i];
        }

	/* free memory from c version */
	free_cubic(grid);
}
