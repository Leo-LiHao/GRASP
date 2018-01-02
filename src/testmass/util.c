/* util.c
   Provides routines not directly associated with waveform generation,
   such as reading files, integrating fuctions etc. 
   Author: S. Droz (droz@physics.uoguelph.ca)
*/
static char rcsid[] = "$Id: util.c,v 1.2 1998/06/06 20:47:23 ballen Exp $\
\n$Name: RELEASE_1_9_8 $";

/* Changes:
   $Log: util.c,v $
   Revision 1.2  1998/06/06 20:47:23  ballen
   Serge modified this code.

   Revision 1.5  1998/04/01 21:02:44  droz
   This version is in GRASP 1.6.5

   Revision 1.4  1998/02/26 22:27:26  droz
   Implemented some of Bruce Allen's suggestions.

   Revision 1.3  1998/02/25 21:40:51  droz
   Integrated this version into the grasp library.

   Revision 1.2  1998/01/29 15:01:30  droz
   All routines seem to perform according to specifications.
   We should replace fopen by GRASP_open.

   Revision 1.1  1998/01/23 18:31:03  droz
   Initial revision

   


*/
 



/* Possible Errors: */

                     
/* We want a general purpos integrator to calculate
   \[ \Phi(x) = \int^x_x_0 dy f(y). \]
   We want to be able to start at x_0 and step through 
   increasing values of x.
   
   To achieve this we slightly change the function odeint
   from the numerical recepies.
   
   To simplify matters we actually just simplify odeint a bit,
   makeing it more userfriendly.

*/
/* CPP symbols: */

#include "grasp.h"

/* NR defines: */

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))


/* Prototypes: */
static float (*FoX)(float x);
/* NR prototypes: */
float *vector(long nl, long nh);
void free_vector(float *v, long nl, long nh);
void rkqs(float y[], float dydx[], int n, float *x,
        float htry, float eps, float yscal[], float *hdid, float *hnext,
        void (*derivs)(float, float [], float []));

static const char TMDATA[] = "GRASP_PARAMETERS";
int integrateODE(float ystart[],  int nvar, float *x1, float x2, float eps, 
               float h1, float hmin, void (*derivs)(float, float [], float []))
/* Integreate a system of first order ODEs from x1 to x2. 
   This routine is adapted from the routine odeint out of
   the book "Numerical Recipes in C", page 721. 
   
   Input: ypstart[]: Initial data. This array is updated 
                     so a next step can be taken by simply changing
		     x2 to a new value.

	  nvar:      The number of variables.
	  *x1:       start value.	
	  x2:        end value
	  eps:	     tolerance
	  h1:        initial stepsize (guess !)
	  hmin:      minimum stepsize.
	  *dervis:   a function deris(x,y,dy); the ODE.
  Errors:
  0:    No error
  5:    Step size too small.
  6:    Too many steps.	  
*/	  	     	
{
	int nstp,i;
	float x,hnext,hdid,h;
	float *yscal,*y,*dydx;
        
        const int MAXSTP = 10000;
        const float TINY = 1.0e-25; /* NR gives 1.0e-30 here, but that seems to be
		                               too small for some applications like e.g.
									   lorenz . */
		
	yscal=vector(1,nvar);
	y=vector(1,nvar);
	dydx=vector(1,nvar);
	x=*x1;
	h=SIGN(h1,x2-*x1);

	for (i=1;i<=nvar;i++) y[i]=ystart[i];
	
	for (nstp=1;nstp<=MAXSTP;nstp++) {
		(*derivs)(x,y,dydx);
		for (i=1;i<=nvar;i++)
			yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
		if ((x+h-x2)*(x+h-*x1) > 0.0) h=x2-x;
		rkqs(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs);
		if ((x-x2)*(x2-*x1) >= 0.0) {
			for (i=1;i<=nvar;i++) ystart[i]=y[i];
			free_vector(dydx,1,nvar);
			free_vector(y,1,nvar);
			free_vector(yscal,1,nvar);
			*x1 = x2;
			return kBhptNoError;
		}
		if (fabs(hnext) <= hmin) 
		{
		  GR_start_error("integrateODE()", rcsid,__FILE__, __LINE__);
      	  GR_report_error("Step size too small.\n");
          GR_end_error();
		  return kBhptStepTooSmall;
        }
		h=hnext;
	}
	GR_start_error("integrateODE()", rcsid,__FILE__, __LINE__);
      	GR_report_error("Too many steps needed.\n");
        GR_end_error();
	return kBhptTooManySteps;
}


static void deriFun(float x, float *y, float *dy)
{
   dy[1] = FoX(x); /* NR array conventions used. */
}

int integrate_function(float vl, float vr, float vo, float (*f)(float ),
                        float **F, int number_of_points)
/* Calculate F(v) = \int_vo^v dx f(x) between vl and vr.
   Return the data in **F anc calculate number_of_points datapoints.
   If **F == NULL allocate the memory.
   Errors:
   0: No error.
   1: -
   2: Not enough memory
*/
{
   int jo,i;
   float dx, x, y;
   const float epsilon = 1.0E-7;
   /* Allocate memory if necessary */
   if ( ! *F )  *F = malloc(number_of_points*sizeof(float));
   if (! *F ) 
   {
      GR_start_error("integrate_function()", rcsid,__FILE__, __LINE__);
      GR_report_error("Not enough memory for y.\n");
      GR_end_error();
      return kBhptOutOfMemory;
   }
   
   FoX = f;
   
   dx = (vr-vl)/(number_of_points-1.0);
   jo = (int)floor((vo-vl)/dx);
   if ( jo < 0 || jo > number_of_points ) return 1;
   /* First integreate to the left: */
   x = vo;
   y = 0;
   /* Get from vo to v_jo */
   
   integrateODE(&y-1, 1, &x, vl + (jo+1)*dx, epsilon, dx/200.0, dx*1.0E-15, &deriFun);
   *(*F + jo+1) = y;
   for (i = jo+2; i < number_of_points; i++)
   {
      integrateODE(&y-1, 1, &x, x+dx, epsilon, dx/200.0, dx*1.0E-15, &deriFun);
      *(*F + i ) = y; 
    
   }
   
   /* Now integreate to the right: */
   x = vo;
   y = 0;
   /* Get from vo to v_jo */
   integrateODE(&y-1, 1, &x, vl + jo*dx, epsilon, dx/200.0, dx*1.0E-15, &deriFun);
   *(*F + jo) = y;
   for (i = jo-1; i >= 0; i--)
   {
      integrateODE(&y-1, 1, &x, x-dx, epsilon, dx/200.0, dx*1.0E-15, &deriFun);
      *(*F + i ) = y;

    }
    
   return 0;
 }  
   
/* The following two routines allow you to read the data files needed
   into memory */   

int read_real_data_file(const char *filename, float **x, float **y, int *number_of_points,
     int ReadX)   
/* Read a data file of the form
  x0 y0
  x1 y1
  ...
  
  Input:
  ------
  *filename: 	Pointer to a string containing the name of the data file.
  **x: 			Pointer to a pointer to an array of floats which will contain
       			the x - values (if ReadX != 0).
       			If *x == NULL allocate the memory.
  **y: 			Pointer to a pointer to an array of floats which will contain
       			the y - values.
       			If *y == NULL allocate the memory.
  *number_of_points: The number of data points to be read in. If = 0 read in all
  				the points found in the file. (This requires two accesses to the file).
  ReadX:		If zero do not read in the x-values.
  
  Output:
  -------
  **x:			x-values.
  **y:			y-values.
  *number_of_points: Number of points read.
  
  					
  Errors returned:
  1 - Could not open file
  2 - out of memory
  3 - Not enough data points found
  
*/
{
 int i = 0;
 FILE *fp;
 float a,b;
 
 /* If *number_of_points == 0 we read through the file to get all
    the points available. 
 */
 
 if ( ! *number_of_points ) 
 {
   if ( 0 == (fp = grasp_open(TMDATA,filename,"r") ))
   {
      GR_start_error("read_real_data_file()", rcsid,__FILE__, __LINE__);
      GR_report_error("Could not open file %s.\n",filename);
      GR_end_error();
      return 1;
   }
   while (fscanf(fp,"%f %f",&a,&b) == 2) (*number_of_points)++;
   fclose(fp);
 }
 
 /* Check if we have to allocate memory ourselfs */
 
 if  (ReadX && ! *x) *x = malloc(*number_of_points*sizeof(float));
 if (! *x ) 
 {
      GR_start_error("read_real_data_file()", rcsid,__FILE__, __LINE__);
      GR_report_error("Not enough memory for x.\n");
      GR_end_error();
      return 2;
   }
 if (! *y) *y = malloc(*number_of_points*sizeof(float));
 if (! *y ) 
 {
      GR_start_error("read_real_data_file()", rcsid,__FILE__, __LINE__);
      GR_report_error("Not enough memory for y.\n");
      GR_end_error();
      return 2;
   }
 
 /* We now should have all the memory anc can start to read the data in */
 
 if ( 0 == (fp = grasp_open(TMDATA,filename,"r") )) 
 {
    GR_start_error("read_real_data_file()", rcsid,__FILE__, __LINE__);
    GR_report_error("Could not open file %s.\n",filename);
    GR_end_error();
    return 1;
 }
 while (i < *number_of_points && fscanf(fp,"%f %f",&a,&b) == 2 ) 
 { 
   *(*y + i ) = b;
   if ( ReadX  )  *(*x + i ) = a;
   i++;
 }
 fclose(fp);
 
  if (i < *number_of_points ) 
  {
      GR_start_error("read_real_data_file()", rcsid,__FILE__, __LINE__);
      GR_report_error("Not enough data points found in file %s.\n",filename);
      GR_end_error();
      return 3;
   }
   else return 0;
}
 
int read_modes(const char *filename, float **x, float **ReA, float **ImA, int *number_of_points,
                int *MaxL, int ReadX)   
/* This routine reads the modes A_{lm} from a data file into memory. The file is
assumed to be of the form 
2 1 
v_0 (Re[A_{21}]_0,Im[A_{21}]_0)
v_1 (Re[A_{21}]_1,Im[A_{21}]_1)
...

2 2
v_0 (Re[A_{22}]_0,Im[A_{22}]_0)
v_1 (Re[A_{22}]_1,Im[A_{22}]_1)
...
3 1
....

Each mode must have all it's components in the file. 

number_of_points referes to the number of points of each individual mode. 
MaxL referes to the maximal number of multipoles to be read (i.e. *MaxL = 3 reads 
the modes 2 1, 2 2, 3 1, 3 2, 3 3).

Input:
  ------
  *filename: 	Pointer to a string containing the name of the data file.
  **x: 			Pointer to a pointer to an array of floats which will contain
       			the x - values (if ReadX != 0).
       			If *x == NULL allocate the memory.
  **ReA: 		Pointer to a pointer to an array of floats which will contain
       			the real part of A_{lm}.
       			If *ReA == NULL allocate the memory.
  **ImA: 		Pointer to a pointer to an array of floats which will contain
       			the imaginary part of A_{lm}.
       			If *ImA == NULL allocate the memory.
  *number_of_points: The number of data points to be read in. If = 0 read in all
  				the points found in the file. (This requires two accesses to the file).
  *MaxL:		The number of l's to read. if zero count and read all of them.
  ReadX:		If zero do not read in the x-values.
  
  Output:
  -------
  **ReA:			x-values.
  **ImA:			y-values.
  *MaxL:			l's read.
  *number_of_points: Number of points read.
Errors:
  1 - Could not open file
  2 - out of memory
  3 - Not enough data points found. The returnrd data is useless. Memory allocated is not freed.
  4 - Inconsitent file format

*/
{
  int NoPo, Lm, Mm,l,m, i,j;
  float a,b,c;
  char line[255];
  FILE *fp;
  
  /* First check if we need to count modes and, if yes do so. */
  if (*MaxL == 0) 
     if ( 0 == ( fp = grasp_open(TMDATA,filename,"r") )) 
     {
        GR_start_error("read_modes()", rcsid,__FILE__, __LINE__);
        GR_report_error("Could not open file %s.\n",filename);
        GR_end_error();
        return 1;
     }   
     else
     {
        /* Read the first mode to get the number of points (assumed to be the same for
           each mode). This part only does minimal error checking. A
           more stringent policy is addopted in the reading part below.  */
       fgets(line,254,fp); 
       if ( sscanf(line,"%d %d",&l,&m) != 2 ) 
       {
         GR_start_error("read_modes()", rcsid,__FILE__, __LINE__);
         GR_report_error("File %s seems corrupt.\n",filename);
         GR_end_error();
         fclose(fp); 
	 return 4; 
       }
       NoPo = 0;
       fgets(line,254,fp); 
       while (sscanf(line,"%f (%f,%f)",&a,&b,&c) == 3 ) 
       { 
         ++NoPo;
         fgets(line,254,fp); 
       }
       if ( ! *number_of_points ) *number_of_points = NoPo;
       /* We now should have the string " 2 2" in line. If not, then we deal
          with a corrupt file */
       Lm = 2; Mm = 1;
       while ( sscanf(line,"%d %d",&l,&m) == 2 )
       {  
         if ( Mm == Lm )  { Lm++;  Mm = 0; } /* Mm = 0 because we increment on the next line */
         Mm++;
         fgets(line,254,fp);
         while (sscanf(line,"%f (%f,%f)",&a,&b,&c) == 3 ) 
	   if ( ! fgets(line,254,fp) ) line[0] = 0; /* We need this because
	  the gnu version of fgets does not change line
	  if no input is read. */
       } 
       fclose(fp);
       if ( Lm != Mm ) 
       {
         GR_start_error("read_read_modes()", rcsid,__FILE__, __LINE__);
         GR_report_error("File %s seems corrupt.\n",filename);
         GR_end_error();
	 return 4; 
       }
       *MaxL = Lm;
     }
  /* We now have the maximum number of Multipole mode to read in. Let us check if we
      we know how many data points. */
  if ( ! *number_of_points ) 
  if ( 0 == (fp = grasp_open(TMDATA,filename,"r") )) return 1;   
     else
     {
       if ( fscanf(fp,"%d %d",&l,&m) != 2 ) 
       {
         GR_start_error("read_read_modes()", rcsid,__FILE__, __LINE__);
         GR_report_error("File %s seems corrupt.\n",filename);
         GR_end_error();
         fclose(fp); 
	 return 4; 
       }
       while (fscanf(fp, "%f (%f,%f)",&a,&b,&c) == 3 )  (*number_of_points)++;
     }  
     
  /* Check if we have to allocate memory */
  if  (ReadX && ! *x) *x = malloc(*number_of_points*sizeof(float));
  if (! *x ) 
  {
      GR_start_error("read_read_modes()", rcsid,__FILE__, __LINE__);
      GR_report_error("Not enough memory for x.\n");
      GR_end_error(); 
      return 2; 
  }
  if (! *ReA) *ReA = malloc(*number_of_points*sizeof(float)* ((*MaxL-1)*(*MaxL+2))/2);
  if (! *ReA ) 
  {
      GR_start_error("read_read_modes()", rcsid,__FILE__, __LINE__);
      GR_report_error("Not enough memory for ReA.\n");
      GR_end_error(); 
      return 2; 
  }
  
  if (! *ImA) *ImA = malloc(*number_of_points*sizeof(float)* ((*MaxL-1)*(*MaxL+2))/2);
  if (! *ImA ) 
  {
      GR_start_error("read_read_modes()", rcsid,__FILE__, __LINE__);
      GR_report_error("Not enough memory for ImA.\n");
      GR_end_error(); 
      return 2; 
  }


  /*  We are now ready to start reading the data into the memory. We 
      do quite a bit of error checking here to be sure we don't read
      any garbagge in. */
  i = 0;
  if ( 0 == (fp = grasp_open(TMDATA,filename,"r") )) return 1; 
  fgets(line,254,fp);
  for (l = 2; l <= *MaxL; l++)  
    for (m = 1; m<= l; m++)
    { 
      if ( sscanf(line,"%d %d",&Lm,&Mm) != 2 || Lm != l || Mm != m) { fclose(fp); return 4; }
      j = 0;
      fgets(line,254,fp);
      while ( sscanf(line,  "%f (%f,%f)",&a,&b,&c) == 3 )
      {
        if (j < *number_of_points)
        {
          if (  ReadX ) 
            if ( i == j) *(*x + j) = a;
            else if (*(*x + j) != a) 
	    {
              GR_start_error("read_read_modes()", rcsid,__FILE__, __LINE__);
              GR_report_error("File %s seems corrupt.\n",filename);
              GR_end_error();
              fclose(fp); 
	      return 4; 
            }
          
          *(*ReA + i) = b;
          *(*ImA + i) = c; 
          i++;
        }
        j++;
        if ( ! fgets(line,254,fp) ) line[0] = 0;
      }
      
      if ( j < *number_of_points) 
      {
         GR_start_error("read_real_data_file()", rcsid,__FILE__, __LINE__);
         GR_report_error("Not enough data points found in file %s.\n",filename);
         GR_end_error();
         fclose(fp); 
	 return 3; 
      }
    }
  fclose(fp);
  return 0;
      
}       


/*********************************************************************
  minustwoYlm
  ********************************************************************/

float minustwoSlm(float theta, int l, int m)
/* Calculate {}_{-2}Y_{lm}(\theta, 0) using the GRASP function
   sw_spheroid. 
*/
{
   float eigenvalues[4] = {0,0,0,0};
   const float norm = 1.0/sqrt(2.0*M_PI);
   float r,i;
    
   eigenvalues[2] = (l + 2)*(l - 2 + 1);
   
   sw_spheroid(&r, &i, cos(theta), 1, 0.0, -2, l, m, eigenvalues);   
   return r*norm;
}
