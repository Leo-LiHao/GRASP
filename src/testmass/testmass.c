/* testmass.c
   Routines to calculate the phase function for 
   the waveforms from black hole perturbation theory.
   Author: S. Droz (droz@physics.uoguelph.ca)
*/
static char rcsid[] = "$Id: testmass.c,v 1.2 1998/06/06 20:47:22 ballen Exp $\
\n$Name: RELEASE_1_9_8 $";

/* Changes:
   $Log: testmass.c,v $
   Revision 1.2  1998/06/06 20:47:22  ballen
   Serge modified this code.

   Revision 1.10  1998/04/01 21:08:49  droz
   Fixed a sign error in the calculation of hcross.

   Revision 1.9  1998/04/01 21:04:08  droz
   This version is in GRASP 1.6.5

   Revision 1.8  1998/02/26 22:26:49  droz
   Implemented some of Bruce Allen's suggestions.

   Revision 1.7  1998/02/25 21:40:18  droz
   Integrated this version into the grasp library.

   Revision 1.6  1998/02/16 20:01:06  droz
   This version seems to pass the obvious tests with different
   memory allocation schemes.

   Revision 1.5  1998/02/16 14:41:57  droz
   This version works, at least if all the memory management is done
   automatically.

   Revision 1.4  1998/02/09 22:30:12  droz
   This version calculates v(t) much more accurately and faster.
   It still needs cleaning up and removing of debbuging code.

   Revision 1.3  1998/02/06 14:31:07  droz
   This version calculates the wavefor (supposedly) correct. However there is
   a problem calculating v(t). Since the DE for v(t) becomes singular at
   v = 1/sqrt(6) there is a numerical problem. This will be fixed in the
   next version.

   This version also allows for arbitraty modes to be in/excluded.

   Revision 1.2  1998/01/29 15:00:18  droz
   All the routines seem to work for the simple test model. We still
   need to do a lot more checking, especially of all the memory
   allocation stuff.

   Revision 1.1  1998/01/23 18:31:36  droz
   Initial revision

*/


#include "grasp.h"


/****************************************************************
  defines
 ****************************************************************/
 
#define kInterpolateOrder 4                      /* Order to use in polynomial interpolation */
#define mode(l,m) (((l)*(l-1) >> 1) + m - 2)   /* calculate the array index for (l,m) */
#define mode2(l,m)  ((l)*(l) + m - 2 - ((m > 0) ? 1 : 0)) /* Don't you like C? */
												 /* Mode for l<= m <= l, m!= 0 */
#define sgn(l)  ( 1 - ( ((l) & 1) << 1 ) )       /* sgn(l) = (-1)^l */
#define FofV(v) pow(v,3)/(2.0*M*M_PI)
#define VofF(f) pow((f)*2.0*M_PI*M, 1.0/3.0)
/* These are from nrutil.h */
static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))
static float minarg1,minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
        (minarg1) : (minarg2))
static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))
	  

/****************************************************************
  Local Prototypes
 ****************************************************************/
/* a few Prototypes to keep lint and friends happy. */


static float fluxP(float v);
static float fluxQ(float x);
static float omega(float x);
static float TofV(float v);
static float VofT(float v);
static float Phase(float v, float *Phi);
/*static void  Vdot( float t, float *v, float *dv);*/
static int GetMoreMemory(float **hp, float **hc, float **f);
static void Alm(float v, int l, int m, float *Re, float *Im);
int calculate_T(float **T, float *v, int NoPo);
float minustwoSlm(float theta, int l, int m);


/* The NR routines we use: */
void polint(float xa[], float ya[], int n, float x, float *y, float *dy);
void hunt(float xx[], unsigned long n, float x, unsigned long *jlo);

/****************************************************************
  Global Variables
 ****************************************************************/
/* We need some globals to share some info across functions. 
   All globals start with a g.*/

static float gDv, *gPtrP, *gPtrv, *gPtrReA, *gPtrImA, *gPtrT;
            /* set to 1 if we allocate memory ourselfs for the respective variables  */
static int	gAv = 0;          /* = 1 if we allocate memory for v ourselfs            */
static int	gNoPoints = 0;    /* Number of data points read in from the data files   */
static int  gMaximumL = 0;    /* Highes l available. We require all m <= gMaximumL   */
static int	gDataRead = 0;    /* = 1 if we have all the data needed                  */
static int	gPhaseCalc = 0;   /* = 1 if the Phase has been calculated                */
static int	gTimeCalc = 0;    /* = 1 if the t(v) has been calculated                 */
const int   kNumberOfFloats = 512;  /* Number of floats allocated in one chunk       */

/****************************************************************
  fluxQ
 ****************************************************************/
static float fluxQ(float x)
/* Calculate (1-6v^2)/(1-3v^2)^(3/2) */
{ 
  float vv = x*x;
  return (1.0 - 6.0*vv)*pow(1.0 - 3.0*vv,-1.5);
} /* End fluxQ */

/****************************************************************
  Alm
 ****************************************************************/
static void Alm(float v, int l, int m, float *Re, float *Im)
/* Returns the value of Alm at v, using a polynomial inter-
   polation.
   Right now we assume equally space data points.
   This simplifies things a bit, and Eric's data files
   satisfy this request. */ 
{
   int k,j;
   float dvalA;
   
   /* Calculate the index to give to the interpolation Routine */  
   k = (int)floor((v-*gPtrv)/gDv);
   j = IMIN(IMAX(k-(kInterpolateOrder-1)/2,1),gNoPoints-kInterpolateOrder)-1;
   /* The -1 makes sure that the Numerical Recepies understand or arrays */
   polint(&gPtrv[j],(gPtrReA+j + gNoPoints*mode(l,m)), kInterpolateOrder, v, Re, &dvalA);
   polint(&gPtrv[j],(gPtrImA+j + gNoPoints*mode(l,m)), kInterpolateOrder, v, Im, &dvalA);
   /* Maybe we should check for the error estimate 
      dvalA and take appropriate action if things are bad? */      
} /* End Alm */

/****************************************************************
  fluxP
 ****************************************************************/
static float fluxP(float v)
/* Returns the value of P at v, using a polynomial inter-
   polation.
   See comments for Alm( ). */ 
{  
   int k,l;
   float valP, dvalP;
   
   /* Calculate the index to give to the interpolation Routine */
   k = (int)floor((v-*gPtrv)/gDv);
   l = IMIN(IMAX(k+1-(kInterpolateOrder-1)/2,1),gNoPoints-kInterpolateOrder);
   polint(&gPtrv[l],&gPtrP[l], kInterpolateOrder, v, &valP, &dvalP);
   return valP;
} /* End fluxP */

/****************************************************************
  Phase
 ****************************************************************/
static float Phase(float v, float *Phi)
/* Returns the value of Phi at v,  using a polynomial inter-
   polation of the data sored in *Phi.
   Right now we assume equally space data points.
   We assume the data to lie on the same grid as P and Alm.  */ 
{  
   int k,l;
   float val, dval;
   
   k = (int)floor((v-*gPtrv)/gDv);
   l = IMIN(IMAX(k-(kInterpolateOrder-1)/2,1),gNoPoints-kInterpolateOrder)-1;
   polint(&gPtrv[l],&Phi[l], kInterpolateOrder, v, &val, &dval);
   return val;
} /* End Phase */

/****************************************************************
  omega
 ****************************************************************/
static float omega(float v)
/* This is the \omega(t(v)) = \dot{\Phi}(t(v))*/
{
  return 5.0*fluxQ(v)/(fluxP(v)*pow(v,6)*32.0);
} /* End omega */
#ifdef No /* This function is not needed in this release */
/****************************************************************
  Vdot
 ****************************************************************/
 static void  Vdot( float t, float *v, float *dv)
/* dv /dt = .... 
   Note that the factor mu/M^2 is NOT included. 
   We use NR notation */
{
  dv[1] =  (pow(v[1],9)*fluxP(v[1])*32.0)/( fluxQ(v[1]) * 5.0);
} /* End Vdot */
#endif
/****************************************************************
  tdot
 ****************************************************************/

static float tdot( float v)
/* dv/dt = .... 
   Note that the factor mu/M^2 is NOT included. 
   We use NR notation. To use in an integration routine. */
{
  return (pow(v,-9)*fluxQ(v)*5.0)/( fluxP(v) * 32.0);
} /* End Tdot */

/****************************************************************
  calculate_T
 ****************************************************************/
 
int calculate_T(float **T, float *v, int NoPo)
/* Calculate the time t as a function of the orbital velocity
   v. We can then use interpolation to get either t(v) or v(t).
   This is faster and more accurate than calculating v(t) on the
   fly while generating the wave forms. */		 
{
   int error = 0;
   float Vmax = FMIN(*(gPtrv + NoPo -1), 1.0/sqrt(6.0)-0.01 );
   
   if (! gDataRead ) return kBhptNoDataRead;
   /* Allocate memorry if needed */
   if (! *T) *T = malloc(NoPo*sizeof(float));
   if (! *T ) 
   {
      GR_start_error("calculate_T()", rcsid,__FILE__, __LINE__);
      GR_report_error("Not enough memory for y.\n");
      GR_end_error();
      return kBhptOutOfMemory;
   }
   /* It seems we get better accuracy if we start from the right. */
   error = integrate_function(*gPtrv, *(gPtrv+NoPo-1),Vmax, &tdot, T, NoPo);
   /* The commented code would calculate the same thing by solviang
      the ODE directly. The obove way seems more practical at the moment.
	  But we might change back ...
   dx = (*(v+NoPo-1) - *v)/(NoPo-1.0);
   y = 0;
   **T = x = *v;
   printf("dx = %f, x = %f NoPo = %d\n",dx,x,NoPo);
   for(i=0; i< NoPo; i++)
   {
     error = integrateODE(&y-1, 1, &x, x+dx, 1.0e-7, dx/200.0, 1.0e-4*dx, &Tdot);
     *(*T + i + 1) = y;
   }
   */
   if (! error ) gTimeCalc = 1;
   return error;
   
} /* End calculate_T */

/****************************************************************
  VofT
 ****************************************************************/
static float VofT( float t )
/* Get v(t) through interpolation. Since the 'x' values are not equally 
   space we have to use the NR routine hunt to get the propoer index. */
{
   static unsigned long k = 0;
   
   int l;
   float valV, dvalV;
   
   hunt(gPtrT-1,gNoPoints, t, &k);
  
   l = IMIN(IMAX(k-1-(kInterpolateOrder-1)/2,1),gNoPoints-kInterpolateOrder)-1;
   
   polint(&gPtrT[l],&gPtrv[l], kInterpolateOrder, t, &valV, &dvalV);
   return valV;
}   /* End VofT */
  
/****************************************************************
  TofV
 ****************************************************************/
static float TofV(float v)
/*  Calculate T(v). Same disclaimers as for fluxP etc apply */ 
{
   int k,l;
   float valT, dvalT;
   
   k = (int)floor((v-*gPtrv)/gDv);
   l = IMIN(IMAX(k-(kInterpolateOrder-1)/2,1),gNoPoints-kInterpolateOrder)-1 ;	 
   polint(&gPtrv[l],&gPtrT[l], kInterpolateOrder, v, &valT, &dvalT);
   return valT;
} /* End fluxT */ 

/****************************************************************
  calculate_testmass_phase
 ****************************************************************/

int calculate_testmass_phase(float fo, float M, float **Phi)
/* Calculate the phase needed for the waveform generation. 
   In this version we calulate the same number of data points
   as the we have available for P. This should not really pose
   a problem, since the phase is a relatively smooth function.
*/		 
{
   int error;
   float Vmax = FMIN(*(gPtrv + gNoPoints-1), 1.0/sqrt(6.0));
   
   if (! gDataRead ) return kBhptNoDataRead;
   if ( (VofF(TSOLAR*fo) > Vmax) || (VofF(TSOLAR*fo) < *gPtrv))
   {
	   GR_start_error("calculate_testmass_phase()", rcsid,__FILE__, __LINE__);
	   GR_report_error("Frequency %f out of range [%e, %e]. \n",
	     fo,pow(*gPtrv,3)/(2.0*M*TSOLAR*M_PI),pow(Vmax,3)/(2.0*M*TSOLAR*M_PI));
	   GR_end_error();
	   return kBhptFOutOfRange;
	}
      
   /* Allocate memorry if needed */
   if (! *Phi ) *Phi = malloc(gNoPoints*sizeof(float));
   if (! *Phi ) 
   {
      GR_start_error("calculate_testmass_phase()", rcsid,__FILE__, __LINE__);
      GR_report_error("Not enough memory for y.\n");
      GR_end_error();
      return kBhptOutOfMemory;
   }
   
   error = integrate_function(*gPtrv, *(gPtrv+gNoPoints-1),
                              VofF(TSOLAR*fo), &omega, Phi, gNoPoints);
   if (! error ) gPhaseCalc = 1;
   return error;
} /* End calculate_testmass_phase */

/****************************************************************
  ReadData
 ****************************************************************/
  
int ReadData(char *filenameP,  char *filenameAlm, float **v,
             int *number_of_points)
/* Read the data files for the modes and the flux into memory.
   if either filename is NULL, use the default names. 
   Read in at most *number_of_points. If *number_of_points == 0
   count them (in the flux file).
   We assume that the flux file and the mode file have the same number
   of points. WE also calculate t(v). */
{
   /* Default filenames: */
   char *fileP = "P.data";
   char *fileA = "Alm.data";
   int error;
   char *filename;
   
   /* First set up the globals so the function fluxP will work. 
      We let the the data reading routines allocate the memory */
   gPtrP   = NULL;
   gPtrReA = NULL;
   gPtrImA = NULL;
   gPtrT   = NULL;
    
   if ( ! *v) gAv = 1; /* We have to allocate memory for v ourselfs */
   
   if ( filenameP ) filename = filenameP; else filename = fileP;
   
   error = read_real_data_file(filename, v, &gPtrP, number_of_points, 1);
   /* If there was an error other than not enough points
      (in which case fewerpoints will have to do) return */
   if ( error && error != kBhptNotEnoughPoints ) return error;
   
   gPtrv = *v;      /* Where is v sstored */
   gDataRead = 1; 
   gDv = (*(gPtrv + *number_of_points-1) - *gPtrv)/(*number_of_points-1.0);
   gNoPoints = *number_of_points;
   
   if ( filenameAlm ) filename = filenameAlm; else filename = fileA;
   error = read_modes(filename, &gPtrv, &gPtrReA, &gPtrImA,  number_of_points, &gMaximumL,0);
   
   /* We now have the data files read into memory. We now have to calculate
      v(t) or rather t(v) */
   if (! error ) error = calculate_T( &gPtrT ,  gPtrv, *number_of_points);
      
   /* Here all errors are fatal. */
   if ( error ) gDataRead = 0;   	     
  
   return error;
} /* ReadData */

/******************************************************************
 GetMoreMemory
 ******************************************************************/
static  int GetMoreMemory(float **hplus, float **hcross, float **fre)
/* Called to (re)allocate memory when calculating the waveform */
 {
    float * flptr;
    static size = 0;
    
    size += kNumberOfFloats;
    
    flptr = realloc(*hcross, sizeof(float)*size);
    if ( flptr ) *hcross = flptr;
    else return kBhptOutOfMemory;
    flptr = realloc(*hplus, sizeof(float)*size);
    if ( flptr ) *hplus = flptr;
    else return kBhptOutOfMemory;
    flptr = realloc(*fre, sizeof(float)*size);
    if ( flptr ) *fre = flptr;
    else return kBhptOutOfMemory;
    
    return 0;
 }
/******************************************************************
 Set_Up_Data
 ******************************************************************/
void Set_Up_Data( float *v, float *P, float *T, float *ReA, float *ImA,
                  int num_of_datapoints  )
/* Let the user suppliy the data files for calculation of the waveform */
{
   if ( v ) gPtrv = gPtrv;
   if ( P ) gPtrv = gPtrP;
   if ( T ) gPtrv = gPtrT;
   if ( ReA ) gPtrv = gPtrReA;
   if ( ImA ) gPtrv = gPtrImA;
   if ( num_of_datapoints ) num_of_datapoints = gNoPoints;
   if ( v && P && ReA && ImA ) gDataRead = 1;
   if ( T ) gTimeCalc = 1;
}

/******************************************************************
 Clean_Up_Memory
 ******************************************************************/
void Clean_Up_Memory( float *Phase )
/* Get rid of all the memory we allocated ourselfs */
{
   free( gPtrP );
   free( gPtrReA );
   free( gPtrImA );
   free( gPtrT );
   if ( Phase )
   {
      free( Phase );
      gPhaseCalc = 0;
   }
   if ( gAv ) free( gPtrv );
   gDataRead = 0;
   gTimeCalc = 0;  
}

 
 
/******************************************************************
 testmass_chirp
 ******************************************************************/
int testmass_chirp(float m1, float m2, float theta, float phi,
               float *Ph, float f_start, float f_end, float *f_started,
               float *f_ended,  float dt, float **hplus, float **hcross, 
	           float **fre, int *number_of_points, int MaxL, int *modes)
/* Calculate the signal from BHPT ... */
{
       int NoPo; 
       float to,v,vl,phase, ReA, ImA;
       float Vmax = FMIN(*(gPtrv + gNoPoints-1), 1.0/sqrt(6.0));
       float M = (m1 + m2)*TSOLAR;
       float mu = (m1*m2/M)*TSOLAR*TSOLAR;
       float beta = M/mu;
       float alpha = mu/(M*M);
       float Vstart = VofF(f_start);
       float Vend   = VofF(f_end);
       int   i,j,k,*m,*l,MaxMod;
       float s,c,*Slmm,*Slmp;
       double arg;
       
       int error = 0;
	       
       /* The data must already be read into  memory */
       if ( ! gDataRead ) return kBhptNoDataRead;
       
       /* and the phase must be calculated */
       if ( ! gPhaseCalc ) return kBhptNoPhase;
       
	   /* and and so does t(v) */
       if ( ! gTimeCalc ) return kBhptNoTime;
	   
       /* Should we allocate memory? */
       if (  *number_of_points ) 
       {
          if (! *hplus ) *hplus = malloc(*number_of_points*sizeof(float));
          if (! *hcross ) *hcross = malloc(*number_of_points*sizeof(float));
          if (! *fre ) *fre = malloc(*number_of_points*sizeof(float));
	  if (! ( *hplus && *hcross && *fre ) ) 
	  {
	    GR_start_error("testmass_chirp()", rcsid,__FILE__, __LINE__);
	    GR_report_error("Not enough memory for the number of points requested..\n");
	    GR_end_error();
	    return kBhptOutOfMemory;
	  }
	} 
	 /* If memory is allocated, but the number of points is not given, 
	    then we are in trouble, since we don't know how much memory we can
	    actually use. */
	if ( ! *number_of_points && ( *hplus || *hcross ) ) 
	{
	   GR_start_error("testmass_chirp()", rcsid,__FILE__, __LINE__);
	   GR_report_error("You need to give the number of points that your memory can hold.\n");
	   GR_end_error();
	   return kBhptUnknownMemory;
	}
	/* A couple of sanity checks */
	*f_started = f_start;
	/* Is the starting frequency large enough ?  */
	if ( Vstart < *gPtrv ) 
	{
	   GR_start_error("testmass_chirp()", rcsid,__FILE__, __LINE__);
	   GR_report_error("Initial %f frequency out of range [%f, %f]. \n",f_start, FofV(*gPtrv),FofV(Vmax));
	   GR_report_error("Adjusted initial frequency\n");
	   GR_end_error();
	   *f_started = FofV(*gPtrv);
	   Vstart = *gPtrv;
	   error = kBhptFOutOfRange;
	}
	if  (Vend > Vmax ) 
	{
	   GR_start_error("testmass_chirp()", rcsid,__FILE__, __LINE__);
	   GR_report_error("Final %f frequency out of range [%f, %f]. \n",f_end, FofV(*gPtrv),FofV(Vmax));
	   GR_report_error("Adjusted final frequency\n");
	   GR_end_error();
	   Vend = Vmax;
	   error =  kBhptFOutOfRange;
	}
	/* Establish which modes will be used. */
	/* We now create an array which containes a table of modes to be used. */
	MaxMod = 0;
	/* How many entries do we need ? */
	for (i = 2; i <= MaxL;i++)
	  for (j=1;j<=i;j++)  if (modes[mode2(i,j)] || modes[mode2(i,-j)] ) MaxMod++; 
	l = malloc(MaxMod*sizeof(int));
	m = malloc(MaxMod*sizeof(int));
	Slmp = malloc(MaxMod*sizeof(float));  
	Slmm = malloc(MaxMod*sizeof(float));
	if ( ! ( l && m && Slmp && Slmm ) ) 
	{ 
	    GR_start_error("testmass_chirp()", rcsid,__FILE__, __LINE__);
	    GR_report_error("Not enough memory for to save mode information.\n");
	    GR_end_error();
	    return kBhptOutOfMemory;
	}
	k = 0;
	for ( i = 2; i <= MaxL; i++) 
	  for (j = 1; j <= i;j++) 
	  if (modes[mode2(i,j)] || modes[mode2(i,-j)] )
	  {
	    l[k] = i; 
	    m[k] = j; 
	    /* There are two combinations we need: Ylm + (-1)^l Yl-m and Ylm - (-1)^l Yl-m.
	       Depending on the modes inluded only the m +- terms should get in */
	    Slmm[k]    = Slmp[k] = 0.0;
        if (modes[mode2(i,j)] )
	      Slmm[k] = -( Slmp[k] = minustwoSlm(theta, i,j)); 
	    if (modes[mode2(i,-j)] )  
	    {
	      Slmp[k] += sgn(i)*minustwoSlm(theta, i,-j);
	      Slmm[k] += sgn(i)*minustwoSlm(theta, i,-j);
	    } 
	    k++;
	  }
	/* Mode stuff finished */ 
	 
	/* Get ready to enter the main loop. */
	/* **fre = f_start; *//* Where we start: v(t=0) = v */
	vl = v = VofF(f_start); 
	to = TofV(vl)/alpha;  

	/*-----------------------------------------------------*/
	/************ Start the main loop **********************/
	
	 for (NoPo = 0 ; v <=Vmax; NoPo++)
	 {
	    if ( *number_of_points )
	    {
	      if (NoPo == *number_of_points ) {NoPo--;break;}
	    } 
	    else 
	    if (  ( NoPo % kNumberOfFloats) ==0  ) /* This we means allocate memory dynamically ! */ 
	      if ( GetMoreMemory(hcross, hplus, fre)  ) 
	      {
	        GR_start_error("testmass_chirp()", rcsid,__FILE__, __LINE__);
	        GR_report_error("Out of Memory.\n");
	        GR_end_error();
     		break;
	      } /* end if */
	 	      
	 phase = Phase(v, Ph);   /* Calculate the Phase $\Phi$ */
	   
	 *(*fre + NoPo ) = FofV(v); /* Convert velocity int frequency */
	 *(*hplus + NoPo) = 0.0;
	 *(*hcross + NoPo) = 0.0;
	 for (i = 0; i< MaxMod; i++)  /* Sum over the modes */
	 {
	    Alm(v, l[i], abs(m[i]), &ReA, &ImA);	
	    arg = m[i] * ( phi - beta * phase); /*since arg can be very large 
	                         it has to be a double, or sin(arg) is shit. */
	    s = sin(arg);
	    c = cos(arg);
	    *(*hplus  + NoPo) += (ReA*c - ImA*s) * Slmp[i];
	    *(*hcross + NoPo) += (ReA*s + ImA*c) * Slmm[i];  
	 }
	 vl = v;
	   
	 /* Get the next v value. */ 
	 v = VofT((to + (NoPo+1) * dt)*alpha); 
     if (v>Vend) break;
   } /* end for */
   *number_of_points = ++NoPo;
   *f_ended = FofV(vl);
   free(m);
   free(l);   
   free(Slmp);
   free(Slmm); 
   return error;
}          
/******************************************************************
 Get_Duration
 ******************************************************************/
float Get_Duration(float f1, float f2, float m1,float m2)
/* Calculate the duration of a chirp betwen tow frequencies */
{
   float v1;
   float v2;
   float M = (m1 + m2)*TSOLAR;
   float mu = (m1*m2/M)*TSOLAR*TSOLAR;
   float alpha = mu/(M*M);

   /* Note that the value M is used in the macro VofF */
   
   if ( ! gTimeCalc ) return -kBhptNoTime;
   
   v1 = VofF(FMIN(f1,f2));
   v2 = VofF(FMAX(f1,f2));
   return (TofV(v2) - TofV(v1))/alpha;
}

/******************************************************************
 Get_Fmax
 ******************************************************************/
float Get_Fmax(float m1,float m2)
/* Calculate maximum (orbital) frequency  */
{
   float M = (m1 + m2)*TSOLAR;
   float Vmax = FMIN(*(gPtrv + gNoPoints-1), 1.0/sqrt(6.0));
   
   /* The data must already be read into  memory */
   if ( ! gDataRead ) return kBhptNoDataRead;   
   /* Note that the value M is used in the macro FofV */
   return (FofV(Vmax));
}
   
