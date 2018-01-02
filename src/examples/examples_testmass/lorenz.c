/* Example lorenz
   Solve the Lorenz equations:
   
   dx/dt = s * (y - x),
   dy/dt = (r * x - y - xz),
   dz/dt = -b*z + xy,
   starting at the point (0.0, 1.0, 1.0 ).
   Try eg. lorenz 11.0 28.0 2.6666 2000 test.dat
   Author: S. Droz
   For more info see for example
   http://pineapple.apmaths.uwo.ca/~blair/lorenzintro.html
   Note that these equations are chaotic, and thus extremely
   susceptible towards numerical errors over long time spans.
   If you change the accuracy (eps below) by a factor of 10 
   you might get a rather different looking picture.
*/ 
#include "grasp.h"

float s,r,b;


void lorenz(float x, float y[], float dy[])
/* Set up the system of ordinary DE. Since we use the NR 
   integrator we use arrays going from 1..3, rather than from 
   0..2. */
{
   dy[1] = s * (y[2] - y[1]);
   dy[2] = (r - y[3]) * y[1] - y[2];
   dy[3] = y[1] * y[2] - b*y[3];
   
}

int main(int argc, char *argv[])
{
  int NoPo,i;
  int error = 0;                   /* No errors yet */
  float y[3]  = { 0, 1.0, 1.0 };   /* Initial point */
  float t = 0.0;                   /* Time always starts at 0 */
  float dt = 0.02;                 /* The time step used */
  FILE *fp;                        
  
  if (argc != 6)  /* Do we have enough command line arguments? */ 
  {
     fprintf(stderr,"Usage: %s s r b number_of_Points filename\n" ,argv[0]);
	 fprintf(stderr,"E.g. %s  11.0 28.0 2.6666 2000 test.dat\n" ,argv[0]);
	 fprintf(stderr,"The numbers s, r and b are parameters of the Lorenz equations\n");
	 fprintf(stderr,"See e.g. http://pineapple.apmaths.uwo.ca/~blair/lorenzintro.html for mor info\n");
	 return 1;
  }
  /* get the command line arguments */
  s = atof(argv[1]);
  r = atof(argv[2]);
  b = atof(argv[3]);
  
  NoPo = atoi(argv[4]);
  fp = fopen(argv[5],"w");
  
  /* Open the output file */
  if (! fp ) {  printf("File error.\n");  exit(1); }
  fprintf(fp, "%10.8e %10.8e %10.8e\n",y[0],y[1],y[2]);
  
  /* Now we start at t=0 and integrate for NoPo-1 time steps: */
  for ( i = 1; i< NoPo; i++)
  {  
     /* So starting at t integrate y to t+dt, using the
	    equations implemented in the function lorenz. 
		We start with an initial step of dt/10 and go
		as low as dt*10^{-10}. We require an accuracy of at least 10^{-6},
		See the NR for a detailed explanation of what all these numbers mean. */
     error = integrateODE(y-1, 3 , &t, t+dt, 1.0e-6, dt/10.0, dt*1.0e-10,
	 &lorenz);
     if (error) break; /* We chicken out if something goes wrong. */
	 fprintf(fp, "%10.8e %10.8e %10.8e\n",y[0],y[1],y[2]);
  }
  fclose(fp);
  return error; /* Good bye */
}  	          
