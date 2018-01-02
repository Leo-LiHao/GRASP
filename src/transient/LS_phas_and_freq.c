/* GRASP: Copyright 1997,1998  Bruce Allen */
/*
 * Integration routine for phase and frequnency functions for Lai and Shapiro
 * waveforms [ApJ 442, 259 (1995)]. Note that we use u := f/fmax for the 
 * frequency variable and that other variables are defined in such a way that 
 * the equations take the form
 *
 * u' = A u^(5.2) [u-1],   
 * 
 * Phi' = 2 pi fmax u
 *
 * where ' denotes differentiation wrt time.
 *
 *
 * Authors:  Warren G Anderson,  warren@ricci.phys.uwm.edu
 *           Patrick R Brady,    patrick@tapir.caltech.edu
 *
 */
#include "grasp.h"

static char *rcsid="$Id: LS_phas_and_freq.c,v 1.3 1998/06/26 19:01:06 ballen Exp $\
\n$Name: RELEASE_1_9_8 $";

/* set initial value for f/fmax */
#define U_INIT 0.9



/* These external variables are used by odeint() */
  int kmax=0,kount=0;
  float *xp=NULL,**yp=NULL,dxsav=0.0;

/* define integration routine */
void LS_phas_and_freq(double Phi[], float u[], float fmax, 
		float dt, int n_samples)
{
  /* declarations */
  int i,nok,nbad;
  float *u0,t0,t1;
  float u0store;

   /* prototypes for Numerical Recipes ODE integrators */
   void odeint(float*, int,float,float,float,float,float,int*,
     int*, void (*derivs)(float,float [],float []), void (*rkqs)(float
     [],float [],int,float*,float,float, float [],float*,float*,void
     (*)(float,float [],float [])));

   void rkqs(float*,float*,int,float*, float,float,float*,
     float*,float*, void (*derivs)(float,float [],float []));

   /* test that memory has been allocated */
  if (Phi==NULL || u==NULL) {
          GR_start_error("LS_phas_and_freq()",rcsid,__FILE__,__LINE__);
          GR_report_error("First two arguments can not be NULL pointers!\n");
          GR_end_error();
          abort();
  }
  
 /* initial values for DE integration */
  Phi[0]=0.0;     /* take initial phase to be 0 */
  u[0]=U_INIT;    /* start integration at 90% of maximum frequency */
  t0=0.0;         /* first integration starts at time 0 */
  t1=dt;          /* first integration ends at time dt */
  u0=&u0store-1 ; /* point to memory for initial values */ 
  u0[1]=U_INIT;   /* assign initial value for first integration */


  /* integrate u equation using adaptive step size 4th order Runge-Kutta 
     and (odeint from numerical recipes). Since Phi is a simple integral
     of u and u changes so slowly, trapezoidal integration is sufficient
     for Phi equation.*/
  for(i=1;i<n_samples;i++){
     odeint(u0,1,t0,t1,1.0e-7,1.0e-3,1.0e-6,&nok,&nbad,LS_freq_deriv,rkqs);
     u[i]=u0[1];
     Phi[i]=Phi[i-1]+M_PI*fmax*(u[i-1]+u[i])*dt;
     t0=t1;
     t1+=dt;
  } 
}


