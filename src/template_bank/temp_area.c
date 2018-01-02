/* GRASP: Copyright 1997,1998  Bruce Allen */
static char *rcsid="$Id: temp_area.c,v 1.2 1998/01/23 17:59:51 ballen Exp $\
\n$Name: RELEASE_1_9_8 $";

#include "grasp.h"
 
/* define a mass ratio limit for the area calculations */
#define MASS_RATIO_LIMIT 1.e+6

/* 
       Area-of-Parameter-Space Calculation
*/
float template_area(struct Scope *Grid) 
{
float area_formula(float massratio);
float y,A1,A2,A3,A30,A0,pfmsun,rmxmn,rmnmx,mx,mn;


/* 

  Routine for finding the area [dimensions: sec^2] of the allowable
  region in the (tau0,tau1) parameter space. The results are obtained
  by explicit analytic integration of equation 3 in Ref. 2.  The
  results are returned to the calling routine as the value of this
  function: template_area(struct Scope *Grid).  I had a real nice
  verion of this routine to find the area based on a numerical
  integration using "qtrap" in Numerical Recipes.

  This is a "stand-alone" routine; it does not require references to
  any numerical-recipe or other functions. It only requires the minimum
  and maximum masses of the system and the cut-off frequency f_0 and
  other defined constants (pi and the solar mass in seconds).
 
  The area is computed by finding the area under the maximum-mass curve
  plus the area under the minimum-mass curve, then subtracting the area
  under the equal-mass curve. Computation of the area under the
  max/min-curves is essentially identical: just carefully interchange
  the roles of maximum and minimum mass.  The area under the equal mass
  curve has a relatively simple formula since tau0 and tau1 can be
  easily written as functions of a single mass parameter running along
  the equal-mass curve.

  The behavior of the two coalescence times (tau0 and tau1) as
  functions of the masses of the objects is quite pathological; they
  depend on inverse powers of the masses. The area inturn depends on
  the log of mass ratio.  Therefore, for extreme mass ratios,  even
  evaluation of the analytic expressions is tricky.  [Numerical
  integration schemes run into similar problems in the corners of the
  parameter space.] This routine has been tested and shown to give
  solid, reliable results for extended range of mass values:  mass
  ratios up to 10^6.  If this mass ratio is exceeded, a warning is
  printed, but the program is not terminated.
 
*/

pfmsun = M_PI*(Grid->f_start)*TSOLAR;
A0=(18575./49545216.)*pow(pfmsun,-14./3.)*TSOLAR*TSOLAR;
A30=(60875./16515072.)*pow(pfmsun,-14./3.)*pow(2,-2./3.)*TSOLAR*TSOLAR;
mx=Grid->m_mx;
mn=Grid->m_mn;
rmxmn = mn/mx;
rmnmx=1./rmxmn;

/* check to see if in extreme mass-ration limit, if so print warning */
if (rmnmx>MASS_RATIO_LIMIT) { 
    GR_start_error("template_area()",rcsid,__FILE__,__LINE__);
    GR_report_error("WARNING! Area of parameter space calculation is untested\n");
    GR_report_error("for mass ratios exceeding %e \n",MASS_RATIO_LIMIT);
    GR_report_error("Area calculations are suspect\n");
    GR_end_error();
}


    /* area under the maximum-mass curve */
    A1=A0*pow(mx,-8./3.)*(area_formula(1.)-area_formula(rmxmn));
    /* area under the minimum-mass curve */
    A2=A0*pow(mn,-8./3.)*(area_formula(rmnmx)-area_formula(1.));
    /* area under the equal mass curve  */
    A3=A30*pow(mn,-8./3.)*(1-pow(rmxmn,8./3.));
    y = A1 + A2 - A3;
    /* printf("%f\t%f\t%f\n",A1,A2,A3); */
return y;
}



/* 
   Evaluation of the anitderivative encountered in finding the area
   under the minimum- and maximum-mass curves.
*/
#define NUM_CON (924./743.)
float area_formula(float u)
{   float area; 
    double dumlog,dumtan,u3,u32,r3;
    u3=pow(1.+u,1./3.); 
    u32=u3*u3;
    r3=sqrt(3.);
    dumtan=(1.+2.*u3)/r3;
    dumlog=(1.+u3+u32)/(1.-2.*u3+u32);

    area=-.5*(3.+2.*(4.+3.*NUM_CON)*u+(5.+9.*NUM_CON)*u*u)/u/u/u32
      +(9.*NUM_CON-1.)*(atan(dumtan)/r3+log(dumlog)/6.);
    return area;
}
