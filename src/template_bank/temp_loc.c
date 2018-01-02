/* GRASP: Copyright 1997,1998  Bruce Allen */
static char *rcsid="$Id: temp_loc.c,v 1.4 1998/03/25 19:48:11 ballen Exp $\
\n$Name: RELEASE_1_9_8 $";

#include "grasp.h"
#include <assert.h>

/* Maximum number of templates that can be generated */
#define N_ABORT 800000

/* define an accuracy for root searching */
#define RTSAFE_ACCURACY 3.e-5

/* which rootfinder should we use.  Set to 1 to use rtsafe() or 0 to use rootfinder() */
#define RTSAFE 0

float x0a[ N_ABORT],x1a[N_ABORT],tau0a[N_ABORT],tau1a[N_ABORT];
double etaa[N_ABORT],mtota[N_ABORT];

/* function that uses Newton-Raphson and bisection to obtain root
   to double precision accuracy
*/
#define ROOTACC 1.e-15
#define MAXITERATIONS 100

#define HOLEFILL 0  /* fill in empty spots along the equal mass line */

/* prototypes */
double rootfinder(void (*fandfprime)(double,double*,double*),double xl,double xh);
int old_m_and_eta(double tau0,double tau1,double *M,double *eta,double Mmin,double Mmax,double pf);

double rootfinder(void (*fandfprime)(double,double*,double*),double xl,double xh)
{

	double xi,fl,fh,fi,dl,dh,di,dxiold,dxi,nrx;
	double alpha,beta,temp;
	int count=0;

	/* first ensure that xl<=xh */
	if (xl>xh) {
		temp=xl;
		xl=xh;
		xh=temp;
	}

	/* now check that they really bisect the root */
	(*fandfprime)(xl,&fl,&dl);
	(*fandfprime)(xh,&fh,&dh);

	/* if we are on a root, return it */
	if (fl==0.0)
		return xl;
	if (fh==0.0)
		return xh;

	if (fl*fh>0.0) {
		GR_start_error("rootfinder()",rcsid,__FILE__,__LINE__);
		GR_report_error("Fatal error -- root not bracketed.\n");
		GR_report_error("f(%f)=%f and f(%f)=%f\n",xl,fl,xh,fh);
		GR_end_error();
		abort();
	}

	/* calculate the midpoint */
	xi=0.5*(xl+xh);
	dxiold=(xh-xl);

	/* start infinite loop, counting */
	while (++count) {

		/* check to see if we have made too many iterations */
		if (count>MAXITERATIONS) {
			GR_start_error("rootfinder()",rcsid,__FILE__,__LINE__);
			GR_report_error("Warning: exceeded %d (MAXITERATIONS) iterations.\n",MAXITERATIONS);
			GR_report_error("The bisection/Newton-Raphson rootfinder is not converging.\n");
			GR_report_error("Returning (possibly inaccurate) root estimate %f (bracketing range is %f to %f).\n",
				xi,xl,xh);
			GR_end_error();
			return xi;
		}

		/* evaluate the function at the intermediate point */
		(*fandfprime)(xi,&fi,&di);
		if (fi==0.0) {
			return xi;
		}

		/* if the bisection boundaries are very close together, return */
		if (fabs((xh-xl)/(xh+xl))<ROOTACC) {
			alpha=fh/(fh-fl);
			beta=1.0-alpha;
			return (alpha*xl+beta*xh);
		}

		/* set bracketing boundaries */
		if (fl*fi>0) {
			fl=fi;
			xl=xi;
		}
		else {
			fh=fi;
			xh=xi;
		}

		/* estimate Newton-Raphson correction */
		if (di!=0.0) {
			dxi=fi/di;
			nrx=xi-dxi;

			/* If Newton-Raphson claims we are really close, then return */
			if (fabs(dxi/xi)<ROOTACC) {
				return xi;
			}

			/* if Newton-Raphson within range, and step small, then use it */
			if (nrx>xl && nrx<xh && fabs(dxi)<0.5*fabs(dxiold)) {
				xi=nrx;
				dxiold=dxi;
			}
			else {
				/* Newton-Raphson failing, so do bisection step */
				xi=0.5*(xl+xh);
				dxiold=(xh-xl);

			}

		}
		else {
			/* bisection step */
			xi=0.5*(xl+xh);
			dxiold=(xh-xl);
		}
	}
	/* this never happens, but keep compiler from complaining */
	return 0.0;
}

/* Routine to create a grid of templates.  This is based on the work
   contained in "Gravitational waves from coalescing binaries:
   detection strategies and monte carlo estimation of parameters" by
   R.  Balasubramanian, B.S. Sathyaprakash, S.V. Dhurandhar, Phys.
   Rev.  D53 (1996) 3033-3055; Erratum in Phys. Rev. D54 (1996) 1860.
   Preprint number: gr-qc/9508011 (Ref 1).  Also useful is "Search
   templates for gravitational waves from inspiraling binaries:  choice
   of template spacing", by Benjamin J. Owen, Phys. Rev. D53 (1996)
   6749-6761, gr-qc/9511032 (Ref 2).


   The INPUT to template() is a structure Grid defined by Scope. The
   routine uses as inputs the following fields of Grid:

   Grid.m_mn:     Smallest allowed mass of either member of binary system

   Grid.m_mx:     Largest allowed mass of either element of binary system

   Grid.theta:    The metric on the manifold of templates may be
		  diagonalized. The angle theta (radians) is the angle
		  measured counter-clockwise from the tau0 direction to
		  the axis whose size is specified by Grid.dp.  The
                  range of theta is -pi/2 to pi/2.

   Grid.dp:       The first diameter of the ellipse.

   Grid.dq:       The second diameter of the ellipse.

   Grid.f_start:  The templates begin at this gw frequency, and chirp up
	          to very high frequencies.

   The OUTPUT of template is:

   Grid.n_mx:     Set to the number of templates required.

   Grid.start:    Pointer to an array of type Template[Grid.n_mx]
*/

#define NCHECK 10

void template_grid( struct Scope *Grid ) 
{
double mtot_equal(double min_tot,double max_tot,double pf,double theta, double xzero);
float dx0,dx1,pf,th_dum,x0ini,mmx,mmn,m1,m2;
int ncolumn;
float x0mx,t0const,t1const,ctheta,stheta,x0,x1;
int nfilt,i;
float tau0,tau1;
float x1center;
float fraction,increase;
double m_save;
struct Template *data;
int failed,attempts;
float offset,extra,mtotal,mred;
double eta_save,eta,mtot;
float x1ini,x0fin,x1fin,distanceA,distanceB,saveme;

mmx = Grid->m_mx;
mmn = Grid->m_mn;
th_dum = Grid->theta;
pf=M_PI*(Grid->f_start);

/* defines tau0 for one solar mass eta=1/4 ref 2 eqn (3.1) */
t0const = (5.0/64.0)*pow(TSOLAR,-5.0/3.0)*pow(pf,-8.0/3.0);

/* defines tau1 for one solar mass for eta=1/4, ref 2 eqn (3.2) */
t1const = (2435.0/8064.0)*pow(pf,-2.0)/TSOLAR;

/* The angle is always measured to the dp axis of the ellipse in the
CCW direction from the tauo axis.  If this angle is not in the range
from -pi/2 to pi/2, it should be shifted into this range. 
*/

if (th_dum>0.5*M_PI || th_dum<-0.5*M_PI) {
	GR_start_error("template_grid()",rcsid,__FILE__,__LINE__);
	GR_report_error("Warning: the angle theta (defining templates) should be in the range -pi/2 to pi/2.\n");
	GR_report_error("The value theta=%f is being shifted into range by adding N x pi (N integer).\n",th_dum);

	th_dum=fmod(th_dum,M_PI);
	if (th_dum<-0.5*M_PI) th_dum+=M_PI;
	if (th_dum> 0.5*M_PI) th_dum-=M_PI;

	GR_report_error("The new value being used is theta=%f\n",th_dum);
	GR_end_error();

	assert(th_dum<=0.5*M_PI && th_dum>=-0.5*M_PI);
}

/* Put ellipses into standardized form, 0<= th_dum <= pi/2 */
if (th_dum>0.0) {
	dx0=Grid->dp;
	dx1=Grid->dq;
}
else {
	th_dum+=0.5*M_PI;
	dx0=Grid->dq;
	dx1=Grid->dp;
}



/* for closest packing, the ellipses need radii increased by (5/4)^1/2 */
fraction=0.5*sqrt(3.0);
increase=sqrt(0.8);
dx0*=increase;
dx1*=increase;

/* You have two possible ways of choosing the x0 axis 
   Try both techniques, and choose the one that gives
   you more closely-spaced points. */
ctheta=cos(th_dum);
stheta=sin(th_dum);

/* work out coordinates of initial point */
x0ini=ctheta*t0const*pow(2.0*mmx,-5.0/3.0)+stheta*t1const/(2.0*mmx);
mtot=2.0*mmx;
x1ini=-stheta*t0const*pow(mtot,-5.0/3.0)+ctheta*t1const/mtot;

/* work out coordinates of final point */
x0fin=x0ini+dx0;
mtot=mtot_equal(1.e-12,1.e12,pf,th_dum,x0fin);

if (mtot>0.0) {
	x1fin=-stheta*t0const*pow(mtot,-5.0/3.0)+ctheta*t1const/mtot;

	/* and the first distance */
	distanceA=dx0*dx0+(x1fin-x1ini)*(x1fin-x1ini);
}
else
	distanceA=1.e30;

ctheta=cos(th_dum-0.5*M_PI);
stheta=sin(th_dum-0.5*M_PI);

/* work out coordinates of initial point */
x0ini=ctheta*t0const*pow(2.0*mmx,-5.0/3.0)+stheta*t1const/(2.0*mmx);
mtot=2.0*mmx;
x1ini=-stheta*t0const*pow(mtot,-5.0/3.0)+ctheta*t1const/mtot;

/* work out coordinates of final point */
x0fin=x0ini+dx1;
mtot=mtot_equal(1.e-12,1.e12,pf,th_dum-0.5*M_PI,x0fin);

if (mtot>0.0) {

	x1fin=-stheta*t0const*pow(mtot,-5.0/3.0)+ctheta*t1const/mtot;

	/* and the second distance */
	distanceB=dx1*dx1+(x1fin-x1ini)*(x1fin-x1ini);
}
else
	distanceB=1.e30;

/* Now decide which way to orient the stepping procedure */
if (distanceA==1.e30 && distanceB==1.e30) {
		GR_start_error("template_grid()",rcsid,__FILE__,__LINE__);
		GR_report_error("Fatal error!\n");
		GR_report_error("Unable to find a direction in which to lay down template bank.\n");
		GR_end_error();
		abort();
}
else if (distanceB<distanceA) {
	th_dum-=0.5*M_PI;
	saveme=dx0;
	dx0=dx1;
	dx1=saveme;
}

/* Angle to rotate the axes by */
ctheta=cos(th_dum);
stheta=sin(th_dum);

/* check that ellipse diameters are positive */
assert(dx0>0.0 && dx1>0.0);

/*    New coordinates:
 *    x0= cos(theta) tau0 + sin(theta) tau1
 *    x1=-sin(theta) tau0 + cos(theta) tau1
 */

/* smallest x0: intersect equal mass & min mass lines */
x0ini=ctheta*t0const*pow(2.0*mmx,-5.0/3.0)+stheta*t1const/(2.0*mmx);

/* largest x0: intersect equal mass & min mass lines */
x0mx= ctheta*t0const*pow(2.0*mmn,-5.0/3.0)+stheta*t1const/(2.0*mmn);

/* set number of columns to zero */
ncolumn=0 ;

extra=0.0;
/* move right on equal mass curve, up to minimum mass curve */
for (x0=x0ini;x0<x0mx;x0+=fraction*dx0+extra) {

	/* find value of total mass along equal mass curve */
	if (x0==x0ini)
		mtot=2.0*mmx;
	else {
		mtot=mtot_equal(2.0*mmn,2.0*mmx,pf,th_dum,x0);
		if (mtot<0.0) {
			GR_start_error("template_grid()",rcsid,__FILE__,__LINE__);
			GR_report_error("Fatal error! Not moving up along equal mass line...\n");
			GR_end_error();
			abort();
		}
	}

	/* find the other coordinate along equal mass curve */
	x1=-stheta*t0const*pow(mtot,-5.0/3.0)+ctheta*t1const/mtot;

	/* store the coordinates x0, x1 */
	x0a[ncolumn]=x0;
	x1a[ncolumn]=x1;


	/* and corresponding coordinates tau0,tau1,eta,totalmass */
	tau0a[ncolumn]=t0const*pow(mtot,-5.0/3.0);
	tau1a[ncolumn]=t1const/mtot;
	etaa[ncolumn]=0.25;
	mtota[ncolumn]=mtot;

	/* offset every other ellipse for optimal packing */
	extra=0.0;
	if (ncolumn%2) {
		failed=1;
		attempts=0;
		offset=1.0;
		while (failed && attempts++<30) {
			/* see if offset point is acceptable */
			extra=(1.0-offset)*dx0*(1.0-fraction*increase);
			x0=x0a[ncolumn]+extra;
			x1=x1a[ncolumn]+offset*0.5*dx1;
			tau0=ctheta*x0-stheta*x1;
			tau1=stheta*x0+ctheta*x1;
			/* if it is, then update arrays */
			if (tau0>0 && tau1>0 && m_and_eta(tau0,tau1,&m_save,&eta_save,mmn,mmx,pf)) {
				x0a[ncolumn]=x0;
				x1a[ncolumn]=x1;
				tau0a[ncolumn]=tau0;
				tau1a[ncolumn]=tau1;
				mtota[ncolumn]=m_save;
				etaa[ncolumn]=eta_save;
				failed=0;
			}
			offset*=0.75;
		}
		if (failed) {
			/* if we have stepped past the final point, retreat slightly! */
			if (x0>=x0mx)
				x0=x0mx*0.99999;
			mtot=mtot_equal(2.0*mmn,2.0*mmx,pf,th_dum,x0);
			if (mtot<0.0) {
				GR_start_error("template_grid()",rcsid,__FILE__,__LINE__);
				GR_report_error("Fatal error! Not moving up along equal mass line...\n");
				GR_end_error();
				abort();
			}

			x0a[ncolumn]=x0;
			x1=-stheta*t0const*pow(mtot,-5.0/3.0)+ctheta*t1const/mtot;
			x1a[ncolumn]=x1;
			tau0a[ncolumn]=t0const*pow(mtot,-5.0/3.0);
			tau1a[ncolumn]=t1const/mtot;
			etaa[ncolumn]=0.25;
			mtota[ncolumn]=mtot;
		}
	}

	/* increment the number of columns */
	if (++ncolumn>N_ABORT) { 
		GR_start_error("template_grid()",rcsid,__FILE__,__LINE__);
		GR_report_error("Fatal error!\n");
		GR_report_error("Number of filters N_ABORT = %d insufficient.\n",N_ABORT);
		GR_report_error("Please increase N_ABORT\n");
		GR_end_error();
		exit(1); 
	}
}

/* set number of filters so far, then loop over each column */
nfilt=ncolumn;
for (i=0;i<=ncolumn;i++) {
	x0=x0a[i];
	x1=x1a[i];
	/* move up column */
	while (1) {
		x1+=dx1;

		/* save coordinates of center of ellipse */
		x1center=x1;

		/* this is the lower point rather than center of ellipse */
		x1-=0.5*dx1;

		/* find values of tau0 and tau1 */
		tau0=ctheta*x0-stheta*x1;
		tau1=stheta*x0+ctheta*x1;

		/* exit loop if outside of quadrant */
		if (tau0<0.0 || tau1<0.0) break;

		/* if lower point not in triangular region, exit loop */
		if (!m_and_eta(tau0,tau1,&mtot,&eta,mmn,mmx,pf)) break;
 		
 		/* find coordinates of the center of filter */
		x1=x1center;
		x0a[nfilt]=x0;
		x1a[nfilt]=x1;
		tau0a[nfilt]=tau0=ctheta*x0-stheta*x1;
		tau1a[nfilt]=tau1=stheta*x0+ctheta*x1;

		/* correct eta and total mass for these points */
		if (!m_and_eta(tau0,tau1,mtota+nfilt,etaa+nfilt,1.e-6,1.e6,pf)) {
			GR_start_error("template_grid()",rcsid,__FILE__,__LINE__);
			GR_report_error("Warning!\n");
			GR_report_error("Filter %d has un-acceptable parameters\n",nfilt);
			GR_end_error();
		}

		/* and increment the number of filters */
		if (++nfilt>N_ABORT) { 
			GR_start_error("template_grid()",rcsid,__FILE__,__LINE__);
			GR_report_error("Fatal error!\n");
			GR_report_error("Number of filters N_ABORT = %d insufficient.\n",N_ABORT);
			GR_report_error("Please increase N_ABORT\n");
			GR_end_error();
			exit(1); 
		}
		
	}
}

#if (HOLEFILL)
/* Now check along equal mass line for any  "holes" */
/* calculate radii of ellipses */
rad0sqr=dx0;
rad0sqr*=0.5;
rad0sqr/=increase;
rad0sqr*=rad0sqr;
rad0sqr=1.0/rad0sqr;

rad1sqr=dx1;
rad1sqr*=0.5;
rad1sqr/=increase;
rad1sqr*=rad1sqr;
rad1sqr=1.0/rad1sqr;

for (i=0;i<ncolumn-2;i+=2){ 
	test4=1;
	for (j=0;j<NCHECK && test4;j++) {
		/* get the masses of two adjacent templates along the equal mass line */
		mtotal1=mtota[i];
		mtotal2=mtota[i+2];
		mass=0.5*(mtotal1+j*(mtotal2-mtotal1)/NCHECK);
		tau_of_mass(mass,mass,pf,&tau0d,&tau1d);
		myx0=ctheta*tau0d+stheta*tau1d;
		myx1=ctheta*tau1d-stheta*tau0d;
		test1=(myx0-x0a[i]  )*(myx0-x0a[i]  )*rad0sqr+(myx1-x1a[i]  )*(myx1-x1a[i]  )*rad1sqr;
		test2=(myx0-x0a[i+2])*(myx0-x0a[i+2])*rad0sqr+(myx1-x1a[i+2])*(myx1-x1a[i+2])*rad1sqr;
		test3=(myx0-x0a[i+1])*(myx0-x0a[i+1])*rad0sqr+(myx1-x1a[i+1])*(myx1-x1a[i+1])*rad1sqr;

		if (test1>1.0 && test2>1.0 && test3>1.0) {
			mtota[nfilt]=2.0*mass;
			etaa[nfilt]=0.25;
			tau0a[nfilt]=tau0d;
			tau1a[nfilt]=tau1d;
			test4=0;
			/* and increment the number of filters */
			if (++nfilt>N_ABORT) { 
				GR_start_error("template_grid()",rcsid,__FILE__,__LINE__);
				GR_report_error("Fatal error!\n");
				GR_report_error("Number of filters N_ABORT = %d insufficient.\n",N_ABORT);
				GR_report_error("Please increase N_ABORT\n");
				GR_end_error();
				exit(1);
			}
		}
	}
}
#endif   /* (HOLEFILL) */

/* store the template information */
data=(struct Template *)malloc(sizeof(struct Template)*nfilt);
if (data==NULL){
	GR_start_error("template_grid()",rcsid,__FILE__,__LINE__);
	GR_report_error("Fatal error!\n");
	GR_report_error("Unable to allocate memory for %d templates\n",nfilt);
	GR_end_error();
	exit(1); 
}

Grid->n_tmplt=nfilt;
Grid->templates=data;
for (i=0;i<nfilt;i++) {
	data[i].num=i;
	data[i].tau0=tau0a[i];
	data[i].tau1=tau1a[i];
	data[i].eta=etaa[i];
	data[i].mtotal=mtotal=mtota[i];
	m1=0.5*mtota[i]*(1.0-sqrt(1.0-4.0*etaa[i]));
	m2=mtota[i]-m1;
	data[i].m1=m1;
	data[i].m2=m2;
	data[i].mred=mred=etaa[i]*mtota[i];
	data[i].f_lo=Grid->f_start;
	data[i].mchirp=pow(mred*mred*mred*mtotal*mtotal,0.2);
}

	GR_start_error("template_grid()",rcsid,__FILE__,__LINE__);
	GR_report_error("template_grid(): created %d templates in %d columns.\n",nfilt,ncolumn);
	GR_end_error();
return;
}


/* given m1 and m2, this function returns tau0, tau1  */
void tau_of_mass(double m1,double m2,double pf,double *tau0,double *tau1) {
	double eta,m;

	m=(m1+m2);
	eta=m1*m2/(m*m);
	m*=TSOLAR;

	*tau0=5.0/(256.0*eta)*pow(m,-5.0/3.0)*pow(pf,-8.0/3.0);
	*tau1=(5.0/(192.0*eta*m))*(743.0/336.0+11.0*eta/4.0)*pow(pf,-2.0);

	return;
}

	

/* This function searches along the equal mass line for a total
   mass with a given value of xzero.  It searches from min_tot to
   max_tot, and returns the value of the total mass.
*/
/* Warning: uses passone,passtwo,passx0 to communicate! */
double passone,passtwo,passx0;
double mtot_equal(double min_tot,double max_tot,double pf,double theta, double xzero) {
	double tau0factor,tau1factor,returnval,f1,f2,fp;
	void zero_f_and_d(double,double*,double*);

	/* along the equal mass (eta=1/4) line */
	tau0factor=(5.0/64.0)*pow(TSOLAR,-5.0/3.0)*pow(pf,-8.0/3.0);
	tau1factor=(2435.0/8064.0)*pow(pf,-2.0)/TSOLAR;
	passone=tau0factor*cos(theta);
	passtwo=tau1factor*sin(theta);

	/* pass the value of x0 which we want to cross */
	passx0=xzero;

	/* check to see if a root exists */
	zero_f_and_d(min_tot,&f1,&fp);
	zero_f_and_d(max_tot,&f2,&fp);

	/* if no bracketed root, return negative number */
	if (f1*f2>0.0) return -1.0;

#if RTSAFE
	returnval=rtsafe(zero_f_and_d,min_tot,max_tot,RTSAFE_ACCURACY);
#else
	returnval=rootfinder(zero_f_and_d,min_tot,max_tot);
#endif


	return returnval;
}

/* Function vanishes at total mass corresponding to given x0 for eta=1/4 */
/* Warning: uses passone,passtwo,passx0 to communicate! */
void zero_f_and_d(double m,double *f,double *fp) {
 
	/* return the function mx */
	*f=passone*pow(m,-5./3.)+passtwo/m-passx0;

	/* and its first derivative wrt m */
	*fp=(-5.0/3.0)*passone*pow(m,-8./3.)-passtwo/(m*m);
	return;
}


/* This subroutine takes as input the values of tau0 and tau1.  It returns
   the value of M (total mass) and eta appropriate to tau0 and tau1.  It
   looks for binary systems where each object has a mass >=Mmin and each
   object has a mass <=Mmax.  If it finds such a system, then it returns
   1, and sets M and eta appropriately.  If no system satisfying these constraints
   is found, the function returns 0 and does not change either M or eta
*/

int old_m_and_eta(double tau0,double tau1,double *M,double *eta,double Mmin,double Mmax,double pf) {
	double mcrit;
	int test1,test2;
	double totmass,m1,m2;
	double myeta;
	double myc2,myc3,pft;
	int metasub(double,double,double*,double*,double,double,double);

	/* solar mass x pi x f0 */
 	pft=pf*TSOLAR;

	/* two possible totalmass roots are above/below mcrit */
	myc2=47552.0*pow(pft,8.0/3.0)*tau0;
	myc3=16128.0*pft*pft*tau1;
	mcrit=pow(0.6*myc3/myc2,1.5);

	/* test which mass ranges to search */
	test1=test2=0;
	if (mcrit<2.0*Mmin || 2.0*Mmax<mcrit) 
		/* search only totalmass range [2 Mmin,2 Mmax]  */
		test1=metasub(tau0,tau1,&totmass,&myeta,2.0*Mmin,2.0*Mmax,pf);
	else {
		/* search in two subintervals */
		test1=metasub(tau0,tau1,&totmass,&myeta,2.0*Mmin,mcrit,pf);
		test2=metasub(tau0,tau1,&totmass,&myeta,mcrit,2.0*Mmax,pf);
	}

	/* should never find two distinct physically allowed solutions */
	if (test1 && test2) {
		GR_start_error("m_and_eta()",rcsid,__FILE__,__LINE__);
		GR_report_error("Fatal error!\n");
		GR_report_error("Found two mass/eta's for one tau0/tau1!\n");
		GR_end_error();
		exit(1);
	}
	
	/* if no solutions found, return 0 */
	if (!test1 && !test2)
		return 0;

	/* calculate the individual masses... */
	m1=0.5*totmass*(1.0-sqrt(1.0-4*myeta));
	m2=totmass-m1;

	/* ...and if they are not within the desired range, return 0 */
	if (m1<Mmin || m2>Mmax)
		return 0;

	/* acceptable in-range solutions - so output values and return 1 */ 
	*M=totmass;
	*eta=myeta;
	return 1;
}

/* Binaries with eta=1/4 can make old_m_and_eta() think they have eta>1/4
   due to numerical error.  Therefore, if old_m_and_eta() chokes, we assume
   eta=1/4 and obtain two values of M from tau0 and tau1.  If the two M's
   agree, we're really at eta=1/4.  If not, we've just violated mmin or mmax.
*/
int m_and_eta(double tau0, double tau1, double *M, double *eta, double mmin, double mmax, double pf) {
  double try1, try2, avg;

  if (old_m_and_eta(tau0, tau1, M, eta, mmin, mmax, pf))
    return 1;
  try1 = pow(12.8*tau0,-.6)*pow(pf,-1.6)/TSOLAR;
  try2 = 0.3019593253968253/tau1/pow(pf,2)/TSOLAR;
  avg = 0.5*(try1+try2);
  if (fabs(0.5*(try1-try2)/avg) > 1.e-6 || avg > mmax || avg < mmin)
    return 0;
  *M = avg;
  *eta = 0.25;
  return 1;
}

/* This routine takes values of tau0 and tau1, and looks to see if
   there is a root combination totalmass and eta where totalmass is in
   the range totm_min to totm_max and eta<=1/4.  If so it returns 1 and sets *M
   and *eta, else it returns zero and does not change *M or *eta.
   This routine sets c1,c2,c3 and in this way communicates with is_root.
*/
#define ETAERROR 1.e-14
static double c1,c2,c3;
int metasub(double tau0,double tau1,double *totmass,double *eta,double totm_min,double totm_max,double pf) {
	double value1,value2,deriv1,deriv2,mymass,pft;
	double myeta;
	void is_root(double,double*,double*);

	/* solar mass x pi x f0 */
 	pft=pf*TSOLAR;

	/* these constants are passed to is_root */
	c1=1155.0*TSOLAR;
	c2=47552.0*pow(pft,8.0/3.0)*tau0;
	c3=16128.0*pft*pft*tau1;

 	/* test values at the two endpoints */
	is_root(totm_min,&value1,&deriv1);
	is_root(totm_max,&value2,&deriv2);

	/* if no root between test values, please return */
	if (value1*value2>0.0)
		return 0;

	/* now find the root accurately */
#if RTSAFE
	mymass=rtsafe(is_root,totm_min,totm_max,RTSAFE_ACCURACY);
#else
	mymass=rootfinder(is_root,totm_min,totm_max);
#endif


	/* calculate eta now that one equation is solved */
	myeta=(5.0/256.0)*pow(mymass*pft,-5.0/3.0)/(pf*tau0);

	/*  return 1 if root is physical */
	if (myeta<=0.25) {
		*totmass=mymass;
		*eta=myeta;
		return 1;
	}
	else if (myeta<=0.25+ETAERROR) {
		*totmass=mymass;
		*eta=0.25;
		return 1;
	}

	else
		return 0;
}

/* The function f vanishes when totalmass is that value appropriate to tau0, tau1 */
/* requires external variables c1,c2,c3 */
void is_root(double totmass,double *f,double *fprime) {

	/* function whose root is totalmass */
	*f=c1+totmass*(c2*pow(totmass,2.0/3.0)-c3);

	/* derivative of that function w.r.t. totalmass */
	*fprime=(5.0/3.0)*c2*pow(totmass,2.0/3.0)-c3;
	return;
}




