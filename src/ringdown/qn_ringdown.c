/* GRASP: Copyright 1997,1998  Bruce Allen */
/* This code written by Jolien Creighton, jolien@tapir.caltech.edu */

#include "grasp.h"
extern double log(double);

static char *rcsid="$Id: qn_ringdown.c,v 1.4 1998/01/23 17:59:48 ballen Exp $\
\n$Name: RELEASE_1_9_8 $";

/* preprocessor macros */
#define CZERO Complex(0,0)
#define CUNIT Complex(1,0)
#define CIMAG Complex(0,1)
#define CRE(x) Complex((x),0)
#define CIM(y) Complex(0,(y))
#define CNEG(z) RCmul(-1,(z))
#define CINV(z) Cdiv(CUNIT,(z))
#define TWOPI (2.0*M_PI)
#define TINY 1.e-30
#define EPS 1.e-7
#define SMALL 1.e-6
#define MAXIT 500
#define NEIGEN 50
#define CTINY Complex(TINY,0)
static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))
static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

/* type definitions */
typedef struct FCOMPLEX {float r,i;} fcomplex;
typedef struct QNMODE {int l,m,n,s,sign; float a;} qnmode;
typedef struct QNEIGEN {fcomplex A,omega;} qneigen;
typedef struct QNEIGENTAB {int l,m,s,num; qneigen *pos,*neg;} qneigentab;

/* complex function declarations */
fcomplex Cprod(int n, ...);
fcomplex Csum(int n, ...);
fcomplex Cexp(fcomplex z);
float Csqmod(fcomplex z);
float Carg(fcomplex z);
/* complex function declarations  */
fcomplex Cadd(fcomplex a, fcomplex b);
fcomplex Csub(fcomplex a, fcomplex b);
fcomplex Cmul(fcomplex a, fcomplex b);
fcomplex Complex(float re, float im);
fcomplex Conjg(fcomplex z);
fcomplex Cdiv(fcomplex a, fcomplex b);
float Cabs(fcomplex z);
fcomplex Csqrt(fcomplex z);
fcomplex RCmul(float x, fcomplex a);
/* qn function declarations */
fcomplex qn_sphwf(float mu, qnmode mode, qneigen eigen);
fcomplex qn_sphwf0(float mu, fcomplex norm, qnmode mode, qneigen eigen);
float qn_sphwf0_sqmod(float mu);
fcomplex qn_sphwf_norm(qnmode mode, qneigen eigen);
fcomplex Ccfrac(void (*coef)(int, fcomplex *, fcomplex *));
void qn_sph_Ccfrac_coef(int n, fcomplex *a, fcomplex *b);
void qn_ang_Ccfrac_coef(int n, fcomplex *a, fcomplex *b);
void qn_rad_Ccfrac_coef(int n, fcomplex *a, fcomplex *b);
void qn_sph_coef(int n, qnmode mode, qneigen eigen,
		 fcomplex *alp, fcomplex *bet, fcomplex *gam);
void qn_ang_coef(int n, qnmode mode, qneigen eigen,
		 fcomplex *alp, fcomplex *bet, fcomplex *gam);
void qn_rad_coef(int n, qnmode mode, qneigen eigen,
		 fcomplex *alp, fcomplex *bet, fcomplex *gam);
void qn_axi_vecfunc(int n, float x[], float fvec[]);
void qn_sph_vecfunc(int n, float x[], float fvec[]);
qneigen qn_eigen_solve(qneigen guess, qnmode mode);
void qn_make_eigentab(qneigentab *eigentab);
void qn_print_eigentab(qneigentab eigentab);
qneigen qn_eigen_seek(qneigentab eigentab, qnmode mode);
qneigen qn_eigen(qnmode mode);
qnmode qn_mode(float a, int s, int l, int m);
/* metric routines */
int qn_set_grid(float dl, struct qnScope grid, int flag);
/* other routines */
void error(char *fmt, ...);
int warning(char *fmt, ...);

/* global variables */
qnmode _mode;
qneigen _eigen;

/* complex function definitions */
fcomplex Cprod(int n, ...)
{
  int i;
  va_list ap;
  fcomplex prod,fact;

  va_start(ap,n);
  prod = CUNIT;
  for (i=0;i<n;i++) {
    fact = va_arg(ap, fcomplex);
    prod = Cmul(prod,fact);
  }
  va_end(ap);
  return prod;
}

fcomplex Csum(int n, ...)
{
  int i;
  va_list ap;
  fcomplex sum,term;

  va_start(ap,n);
  sum = CZERO;
  for (i=0;i<n;i++) {
    term = va_arg(ap, fcomplex);
    sum = Cadd(sum,term);
  }
  va_end(ap);
  return sum;
}

fcomplex Cexp(fcomplex z)
{
  return RCmul(exp(z.r),Complex(cos(z.i),sin(z.i)));
}

float Csqmod(fcomplex z)
{
  return z.r*z.r + z.i*z.i;
}

float Carg(fcomplex z)
{
  return atan2(z.i,z.r);
}

/* complex functions */

fcomplex Cadd(fcomplex a, fcomplex b)
{
        fcomplex c;
        c.r=a.r+b.r;
        c.i=a.i+b.i;
        return c;
}

fcomplex Csub(fcomplex a, fcomplex b)
{
        fcomplex c;
        c.r=a.r-b.r;
        c.i=a.i-b.i;
        return c;
}

fcomplex Cmul(fcomplex a, fcomplex b)
{
        fcomplex c;
        c.r=a.r*b.r-a.i*b.i;
        c.i=a.i*b.r+a.r*b.i;
        return c;
}

fcomplex Complex(float re, float im)
{
        fcomplex c;
        c.r=re;
        c.i=im;
        return c;
}

fcomplex Conjg(fcomplex z)
{
        fcomplex c;
        c.r=z.r;
        c.i = -z.i;
        return c;
}

fcomplex Cdiv(fcomplex a, fcomplex b)
{
        fcomplex c;
        float r,den;
        if (fabs(b.r) >= fabs(b.i)) {
                r=b.i/b.r;
                den=b.r+r*b.i;
                c.r=(a.r+r*a.i)/den;
                c.i=(a.i-r*a.r)/den;
        } else {
                r=b.r/b.i;
                den=b.i+r*b.r;
                c.r=(a.r*r+a.i)/den;
                c.i=(a.i*r-a.r)/den;
        }
        return c;
}

float Cabs(fcomplex z)
{
        float x,y,ans,temp;
        x=fabs(z.r);
        y=fabs(z.i);
        if (x == 0.0)
                ans=y;
        else if (y == 0.0)
                ans=x;
        else if (x > y) {
                temp=y/x;
                ans=x*sqrt(1.0+temp*temp);
        } else {
                temp=x/y;
                ans=y*sqrt(1.0+temp*temp);
        }
        return ans;
}

fcomplex Csqrt(fcomplex z)
{
        fcomplex c;
        float x,y,w,r;
        if ((z.r == 0.0) && (z.i == 0.0)) {
                c.r=0.0;
                c.i=0.0;
                return c;
        } else {
                x=fabs(z.r);
                y=fabs(z.i);
                if (x >= y) {
                        r=y/x;
                        w=sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
                } else {
                        r=x/y;
                        w=sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
                }
                if (z.r >= 0.0) {
                        c.r=w;
                        c.i=z.i/(2.0*w);
                } else {
                        c.i=(z.i >= 0) ? w : -w;
                        c.r=z.i/(2.0*c.i);
                }
                return c;
        }
}

fcomplex RCmul(float x, fcomplex a)
{
        fcomplex c;
        c.r=x*a.r;
        c.i=x*a.i;
        return c;
}

/* qn functions */

fcomplex qn_sphwf(float mu, qnmode mode, qneigen eigen)
{
  fcomplex norm;
  norm = qn_sphwf_norm(mode,eigen);
  return qn_sphwf0(mu,norm,mode,eigen);
}

fcomplex qn_sphwf0(float mu, fcomplex norm, qnmode mode, qneigen eigen)
{
  int n=0,max=MAXIT,s=mode.s,m=mode.m;
  float k1=0.5*abs(m-s),k2=0.5*abs(m+s);
  float fact=1+mu,prod=fact;
  fcomplex sum,alp,bet,gam,a,aa,a_,Delta;

  sum = a_ = norm;
  qn_ang_coef(n,mode,eigen,&alp,&bet,&gam);
  a = Cdiv(Cmul(bet,a_),CNEG(alp));
  while (n++<max) {
    Delta = RCmul(prod,a);
    sum = Cadd(sum,Delta);
    if (Csqmod(Delta)<SMALL) break;
    qn_ang_coef(n,mode,eigen,&alp,&bet,&gam);
    a = Cdiv(Cadd(Cmul(bet,aa=a),Cmul(gam,a_)),CNEG(alp));
    a_ = aa;
    prod *= fact;
  }
  if (n>=max) warning("qn_sphwf0(): maximum iterations=%d exceeded. mu=%f",max,mu);
  sum = RCmul(pow(1+mu,k1)*pow(1-mu,k2),sum);
  sum = Cmul(Cexp(RCmul(mu*mode.a,eigen.omega)),sum);
  return sum;
}

float qn_sphwf0_sqmod(float mu)
{
  return Csqmod(qn_sphwf0(mu,CUNIT,_mode,_eigen));
}

fcomplex qn_sphwf_norm(qnmode mode, qneigen eigen)
{
  float qromb(float (*func)(float), float, float);
  fcomplex norm,S0;
  int sign,q;

  _mode = mode;
  _eigen = eigen;

  S0 = qn_sphwf0(0.0,CUNIT,mode,eigen);
  norm = Cexp(CIM(-Carg(S0)));
  S0 = qn_sphwf0(-1.0+EPS,norm,mode,eigen);
  (S0.r < 0.0) ? (sign = -1) : (sign = 1);
  (mode.m > mode.s) ? (q = mode.m) : (q = mode.s);
  (((mode.l - q)%2)==0) ? (sign *= 1) : (sign *= -1);
  norm = RCmul(sign/sqrt(qromb(qn_sphwf0_sqmod,-1.0,1.0)),norm);
  return norm;
}

fcomplex Ccfrac(void (*coef)(int, fcomplex *, fcomplex *))
{
  int n=0,max=MAXIT;
  fcomplex a,b,f,f_,C,C_,D,D_,Delta;

  (*coef)(n,&a,&b);
  if (Cabs(f=b)<TINY) f = CTINY;
  C = f;
  D = CZERO;
  while(n++<max) {
    f_ = f;
    C_ = C;
    D_ = D;
    (*coef)(n,&a,&b);
    D = Cadd(b,Cmul(a,D_));
    if (Cabs(D)<TINY) D = CTINY;
    C = Cadd(b,Cdiv(a,C_));
    if (Cabs(C)<TINY) C = CTINY;
    D = CINV(D);
    Delta = Cmul(C,D);
    f = Cmul(f_,Delta);
    /* if (Csqmod(Csub(Delta,CUNIT))<SMALL) return f; */
    if (Cabs(Csub(f,f_))<SMALL) return f;
  }
  warning("Ccfrac(): maximum iterations %d exceeded.",max);
  return f;
}

void qn_sph_Ccfrac_coef(int n, fcomplex *a, fcomplex *b)
{
  static fcomplex alp_;
  fcomplex alp,gam;

  if (!n) alp_ = CZERO;
  qn_sph_coef(n,_mode,_eigen,&alp,b,&gam);
  *a = CNEG(Cmul(alp_,gam));
  alp_ = alp;
}

void qn_ang_Ccfrac_coef(int n, fcomplex *a, fcomplex *b)
{
  static fcomplex alp_;
  fcomplex alp,gam;

  if (!n) alp_ = CZERO;
  qn_ang_coef(n,_mode,_eigen,&alp,b,&gam);
  *a = CNEG(Cmul(alp_,gam));
  alp_ = alp;
}

void qn_rad_Ccfrac_coef(int n, fcomplex *a, fcomplex *b)
{
  static fcomplex alp_;
  fcomplex alp,gam;

  if (!n) alp_ = CZERO;
  qn_rad_coef(n,_mode,_eigen,&alp,b,&gam);
  *a = CNEG(Cmul(alp_,gam));
  alp_ = alp;
}

void qn_sph_coef(int n, qnmode mode, qneigen eigen,
		 fcomplex *alp, fcomplex *bet, fcomplex *gam)
{
  fcomplex rho,A=eigen.A,omega=eigen.omega;
  int s=mode.s;

  rho = CNEG(Cmul(CIMAG,omega));
  *alp = Csum(3,RCmul(n,Cadd(CRE(n),RCmul(2,Cadd(CUNIT,rho)))),
	      RCmul(2,rho),CUNIT);
  *bet = CNEG(Csum(4,RCmul(2*n,Csum(3,CRE(n),RCmul(4,rho),CUNIT)),
		   Cmul(RCmul(4,rho),Cadd(RCmul(2,rho),CUNIT)),A,CRE(1+s)));
  *gam = Csum(3,RCmul(n,Cadd(CRE(n),RCmul(4,rho))),RCmul(4,Cmul(rho,rho)),
	      CRE(-s*s));
}

void qn_ang_coef(int n, qnmode mode, qneigen eigen,
		 fcomplex *alp, fcomplex *bet, fcomplex *gam)
{
  fcomplex w,A=eigen.A,omega=eigen.omega;
  int s=mode.s,m=mode.m;
  float a=mode.a,k1=0.5*abs(m-s),k2=0.5*abs(m+s);

  w = RCmul(a,omega);
  *alp = CRE(-2*(n+1)*(n+2*k1+1));
  *bet = Csum(4,CRE(n*(n+2*k1+2*k2+1)+(k1+k2)*(k1+k2+1)-s*(s+1)),
	      RCmul(-2*(2*n+2*k1+s+1),w),CNEG(Cmul(w,w)),CNEG(A));
  *gam = RCmul(2*(n+k1+k2+s),w);
}

void qn_rad_coef(int n, qnmode mode, qneigen eigen,
		 fcomplex *alp, fcomplex *bet, fcomplex *gam)
{
  fcomplex A=eigen.A,omega=eigen.omega;
  float a=mode.a,b=sqrt(1-4*a*a);
  int s=mode.s,m=mode.m;
  fcomplex c0,c1,c2,c3,c4,cfact,rho,omega2;

  rho = CNEG(Cmul(CIMAG,omega));
  omega2 = Cmul(omega,omega);

  cfact = RCmul(1/b,Csub(omega,CRE(2*a*m)));
  c0 = Csum(3,CRE(1-s),rho,CNEG(Cmul(CIMAG,cfact)));
  c1 = Csum(3,CRE(-4),RCmul(-2*(2+b),rho),Cmul(CIM(2),cfact));
  c2 = Csum(3,CRE(s+3),RCmul(3,rho),CNEG(Cmul(CIMAG,cfact)));
  c3 = Csum(6,RCmul(4+2*b-a*a,omega2),RCmul(-2*m*a,omega),CRE(-s-1),
	    RCmul(-2-b,rho),CNEG(A),Cmul(Cadd(RCmul(2,omega),CIMAG),cfact));
  c4 = Csum(4,CRE(s+1),RCmul(-2,omega2),RCmul(2*s+3,rho),
	    CNEG(Cmul(Cadd(RCmul(2,omega),CIMAG),cfact)));
  *alp = Cadd(RCmul(n,Csum(3,CRE(n),c0,CUNIT)),c0);
  *bet = Cadd(RCmul(n,Csum(3,CRE(-2*n),c1,CRE(2))),c3);
  *gam = Csum(4,RCmul(n,Csum(3,CRE(n),c2,CRE(-3))),c4,CNEG(c2),CRE(2));
}

void qn_axi_vecfunc(int n, float x[], float fvec[])
{
  fcomplex cf1,cf2;
  _eigen.omega = Complex(x[1],x[2]);
  _eigen.A = Complex(x[3],x[4]);
  cf1 = Ccfrac(qn_rad_Ccfrac_coef);
  cf2 = Ccfrac(qn_ang_Ccfrac_coef);
  fvec[1] = cf1.r;
  fvec[2] = cf1.i;
  fvec[3] = cf2.r;
  fvec[4] = cf2.i;
}

void qn_sph_vecfunc(int n, float x[], float fvec[])
{
  fcomplex cf;
  int l=_mode.l,s=_mode.s;
  _eigen.omega = Complex(x[1],x[2]);
  _eigen.A = CRE(l*(l+1)-s*(s+1));
  cf = Ccfrac(qn_sph_Ccfrac_coef);
  fvec[1] = cf.r;
  fvec[2] = cf.i;
}

qneigen qn_eigen_solve(qneigen guess, qnmode mode)
{
  void newt(float [], int, int *, void (*vecfunc)(int, float[], float[]));
  qneigen eigen;
  float x[4];
  int check;

  /* assign global variables to pass to root finding routine */
  _mode = mode;
  _mode.m = abs(mode.m);

  /* put initial guess (for m>0) into array */
  x[0] = guess.omega.r;
  x[1] = guess.omega.i;
  x[2] = guess.A.r;
  x[3] = guess.A.i;

  /* solve for the root */
  mode.a == 0.0 ? newt(x-1,2,&check,qn_sph_vecfunc)
                : newt(x-1,4,&check,qn_axi_vecfunc);
  if (check) warning("qn_eigen_solve(): local minimum found.");

  /* store array values in eigen */
  if (mode.m<0) {
    eigen.omega = Complex(-x[0],x[1]);
    eigen.A = Complex(x[2],-x[3]);
  } else {
    eigen.omega = Complex(x[0],x[1]);
    eigen.A = Complex(x[2],x[3]);
  }

  return eigen;
}

void qn_make_eigentab(qneigentab *eigentab)
{
  qnmode mode;
  float fact=1.0/sqrt(27),da=0.5/(*eigentab).num;
  int i,l=(*eigentab).l,m=(*eigentab).m,s=(*eigentab).s,num=(*eigentab).num;

  /* allocate memory if necessary */
  if ((*eigentab).pos==NULL)
    if (((*eigentab).pos=(qneigen *)malloc(num*sizeof(qneigen)))==NULL)
      error("qn_make_eigentab(): could not allocate %d bytes of memory.",num*sizeof(qneigen));
  if ((*eigentab).neg==NULL)
    if (((*eigentab).neg=(qneigen *)malloc(num*sizeof(qneigen)))==NULL)
      error("qn_make_eigentab(): could not allocate %d bytes of memory.",num*sizeof(qneigen));

  /* guesses for the eigenvalues for a=0 */
  (*eigentab).pos[0].omega = RCmul(fact,Complex(2*l+1,-1));
  (*eigentab).neg[0].omega = RCmul(fact,Complex(-2*l+1,-1));
  (*eigentab).pos[0].A = (*eigentab).neg[0].A = Complex(l*(l+1)-s*(s+1),0);

  /* assign the mode fields */
  mode.s = s;
  mode.l = l;
  mode.m = m;
  mode.n = 1;  /* n=1: fundamental */
  mode.a = 0;

  /* construct the table */
  (*eigentab).pos[0] = qn_eigen_solve((*eigentab).pos[0],mode);
  (*eigentab).neg[0] = qn_eigen_solve((*eigentab).neg[0],mode);
  for (i=1;i<num;i++) {
    mode.a += da;
    (*eigentab).pos[i] = qn_eigen_solve((*eigentab).pos[i-1],mode);
    (*eigentab).neg[i] = qn_eigen_solve((*eigentab).neg[i-1],mode);
  }

  return;
}

void qn_print_eigentab(qneigentab eigentab)
{
  float a,da=0.5/eigentab.num;
  int i,num=eigentab.num;

  for (i=0,a=0;i<num;i++,a+=da)
    printf("%f\t%f\t%f\t%f\t%f\n",a,
	   eigentab.pos[i].A.r,eigentab.pos[i].A.i,
	   eigentab.pos[i].omega.r,eigentab.pos[i].omega.i);
  for (i=0,a=0;i<num;i++,a+=da)
    printf("%f\t%f\t%f\t%f\t%f\n",a,
	   eigentab.neg[i].A.r,eigentab.neg[i].A.i,
	   eigentab.neg[i].omega.r,eigentab.neg[i].omega.i);
  return;
}

qneigen qn_eigen_seek(qneigentab eigentab, qnmode mode)
{
  qneigen guess;
  int i,num=eigentab.num;

  if ((mode.s != eigentab.s) || (mode.l != eigentab.l) ||
      (abs(mode.m) != eigentab.m) || (mode.n != 1))
    warning("qn_eigenvalues(): mode does not agree with eigenvalue table.");

  i = IMIN((int)(2*mode.a*num),num-1);
  (mode.sign > 0) ? (guess = eigentab.pos[i]) : (guess = eigentab.neg[i]);
  return qn_eigen_solve(guess,mode);
}

qneigen qn_eigen(qnmode mode)
{
  static qneigentab eigentab;
  static int first=1;

  if (first) { /* create eigen table */
    first = 0;
    eigentab.s = mode.s;
    eigentab.l = mode.l;
    eigentab.m = abs(mode.m);
    eigentab.num = NEIGEN;
    eigentab.pos = NULL;
    eigentab.neg = NULL;
    qn_make_eigentab(&eigentab);
  }

  if ((mode.s != eigentab.s) || (mode.l != eigentab.l) ||
      (abs(mode.m) != eigentab.m) || (mode.n != 1)) { /* reset eigen table */
    eigentab.s = mode.s;
    eigentab.l = mode.l;
    eigentab.m = abs(mode.m);
    qn_make_eigentab(&eigentab);
  }

  return qn_eigen_seek(eigentab,mode);
}

qnmode qn_mode(float a, int s, int l, int m)
{
  qnmode mode;
  mode.l = l;
  mode.m = m;
  mode.n = 1; /* fundamental mode */
  mode.s = s;
  (a<0) ? (mode.sign=-1) : (mode.sign=1);
  mode.a = 0.5*fabs(a); /* factor of 0.5 to convert to Leaver */
  return mode;
}

void qn_eigenvalues(float eigenvalues[], float a, int s, int l, int m)
{
  qneigen eigen;
  qnmode mode;

  mode = qn_mode(a,s,l,m);
  eigen = qn_eigen(mode);
  eigenvalues[0] = 0.5*eigen.omega.r;/* factor of 0.5 to convert from Leaver */
  eigenvalues[1] = 0.5*eigen.omega.i;/* factor of 0.5 to convert from Leaver */
  eigenvalues[2] = eigen.A.r;
  eigenvalues[3] = eigen.A.i;
}

void sw_spheroid(float *re, float *im, float mu, int reset,
		 float a, int s, int l, int m, float eigenvalues[])
{
  static fcomplex norm;
  static qneigen eigen;
  static qnmode mode;
  fcomplex spher;

  if (reset) {
    mode = qn_mode(a,s,l,m);
    /* Leaver eigenvalues: factor of two to convert to Leaver's convention */
    eigen.omega = RCmul(2,Complex(eigenvalues[0],eigenvalues[1]));
    eigen.A = Complex(eigenvalues[2],eigenvalues[3]);
    norm = qn_sphwf_norm(mode,eigen);
  }

  spher = qn_sphwf0(mu,norm,mode,eigen);
  *re = spher.r;
  *im = spher.i;
  return;
} 

int qn_ring(float theta, float phi,
            float eps, float M, float a, int l, int m,
            float du, float atten, int max,
	    float **plusPtr, float **crossPtr)
{
  qnmode mode;
  qneigen eigen;
  fcomplex omega,strain,S,c1,c2;
  int i,ndat=0,s=-2;
  float u,mu=cos(theta),efold=0.1*atten*log(10);

  if (l<abs(s)) error("qn_ring(): l=%d < s=%d",l,abs(s));
  if (abs(m)>l) error("qn_ring(): abs(m=%d) > l=%d",m,l);
  if (fabs(a)>=1) error("qn_ring(): abs(a=%f) >= 1",a);
  if ((*plusPtr)==NULL)
    if (((*plusPtr)=(float *)malloc((size_t)(max*sizeof(float))))==NULL)
      error("qn_ring(): could not allocate %d bytes of memory.",max*sizeof(float));
  if ((*crossPtr)==NULL)
    if (((*crossPtr)=(float *)malloc((size_t)(max*sizeof(float))))==NULL)
      error("qn_ring(): could not allocate %d bytes of memory.",max*sizeof(float));

  mode = qn_mode(a,s,l,m);
  eigen = qn_eigen(mode);
  omega = RCmul(0.5,eigen.omega);

  S = qn_sphwf(mu,mode,eigen);
  c1 = RCmul(-4.0*sqrt(-omega.i*eps)*M/Cabs(omega),Cmul(S,Cexp(CIM(m*phi))));
  c2 = RCmul(-1.0/(M*TSOLAR),Cmul(CIMAG,omega));
  ndat = IMIN(max,1+(int)(efold*M*TSOLAR/(-omega.i*du)));
  for (i=0,u=0;i<ndat;i++,u+=du) {
    strain = Cmul(c1,Cexp(RCmul(u,c2)));
    (*plusPtr)[i] = strain.r;
    (*crossPtr)[i] = -strain.i;
  }
  return ndat;
}

int qn_qring(float psi0, float eps, float M, float a,
             float du, float atten, int max, float **strainPtr)
{
  int i,n;
  float f,g,A,beta,omega,u,efold=0.1*atten*log(10);

  if ((*strainPtr)==NULL)
    if (((*strainPtr)=(float *)malloc((size_t)(max*sizeof(float))))==NULL)
      error("qn_qring(): could not allocate %d bytes of memory.",max*sizeof(float));

  f = 1 - 0.63*pow(1-a,0.3);
  g = pow(1-a,0.45);
  A = -8*M*sqrt(eps*g/(5*f*(16+g)));  /* factor of 5 due to angle averaging */
  beta = f*g/(4*M*TSOLAR);  /* TSOLAR is one solar mass in seconds */
  omega = f/(M*TSOLAR);

  n = IMIN(max,1+(int)(efold/(du*beta)));
  for (i=0,u=0;i<n;i++,u+=du)
    (*strainPtr)[i] = A*exp(-beta*u)*cos(omega*u+psi0);
  return n;
}

int qn_filter(float freq, float qual,
              float du, float atten, int max, float **filterPtr)
{
  int i,n;
  float u,omega=TWOPI*freq,beta=M_PI*freq/qual,efold=0.1*atten*log(10);

  if ((*filterPtr)==NULL)
    if (((*filterPtr)=(float *)malloc((size_t)(max*sizeof(float))))==NULL)
      error("qn_filter(): could not allocate %d bytes of memory",max*sizeof(float));

  n = IMIN(max,1+(int)(efold/(du*beta)));
  for (i=0,u=0;i<n;i++,u+=du)
    (*filterPtr)[i] = exp(-beta*u)*cos(omega*u);
  return n;
}

void qn_normalize(float *u, float *q, float *r, int n, float *norm)
{
  int i,ii,ir,j,m=n>>1;
  float sum;

  sum = q[0]*q[0]*r[0];
  for (j=1;j<m;j++) {
    ii = 1 + (ir = j << 1);
    sum += (q[ir]*q[ir] + q[ii]*q[ii])*r[j];
  }
  sum += q[1]*q[1]*r[m];
  *norm = 1/sqrt(0.5*sum);
  for (i=0;i<n;i++) u[i] = (*norm)*q[i];
  return;
}

void find_ring(float *h, float *u, float *r, float *o,
               int n, int len, int safe, int *off,
               float *snr, float *mean, float *var)
{
  int i,num=n-len-safe;
  float sum1=0,sum2=0,dx;

  *snr = *mean = 0;
  correlate(o,h,u,r,n);
  /* search data from safe to num=n-len-safe */
  for (i=safe;i<num;i++) {  /* first pass: determine maximum and mean */
    if (*snr<fabs(o[i])) *snr = fabs(o[*off=i]);
    *mean += o[i];
  }
  *mean /= (float)(num);
  for (i=safe;i<num;i++) {  /* second pass: compute variance */
    sum1 += (dx = o[i] - *mean);
    sum2 += dx*dx;
  }
  *var = (sum2 - sum1*sum1/(float)(num))/(float)(num-1);
  return;
}

void qn_inject(float *strain, float *signal, float *response, float *work,
               float invMpc, int off, int n, int len)
{
  void realft(float *, unsigned long, int);

  int i,imax,imin;

  if (invMpc==0.0) return;    /* return if the source is infinitely far */
  invMpc *= 2.0/(float)(n);   /* normalization including inverse FFT */
  imin = IMAX(0,off);         /* begining of the ring */
  imax = IMIN(n,off+len);     /* end of the ring or of the data set */
  if (imax<imin) return;      /* return if the ring is not in the data set */

  /* work is the ring waveform with correct normalization and offset */
  for (i=0;i<imin;i++) work[i] = 0;
  for (i=imin;i<imax;i++) work[i] = invMpc*signal[i-off];
  for (i=imax;i<n;i++) work[i] = 0;

  /* put work into the frequency domain and apply ``whitening'' filter */
  realft(work-1,n,1);
  ratio(work+2,work+2,response+2,n/2-1);
  work[0] = work[1] = 0;  /* set DC and Nyquist frequencies to zero */

  /* return the work to the time domain and add it to the strain */
  realft(work-1,n,-1);
  for (i=0;i<n;i++) strain[i] += work[i];
  return;
}

/* metric routines */

int qn_set_grid(float dl, struct qnScope grid, int flag)
{
  float x,dx,x0,x1,y,dy,y0,y1,ds_eff=sqrt(8.0)*dl;
  float Q,Q2,f;
  int count=0;

  x0 = grid.qual_min;
  x1 = grid.qual_max;
  y0 = log(grid.freq_min);
  y1 = log(grid.freq_max);
  for (x=x0;x<x1;x+=dx) {
    Q2 = SQR(Q = x);
    dx = ds_eff*Q*(1+4*Q2)/sqrt(3+16*SQR(Q2));
    for (y=y0;y<y1;y+=dy) {
      f = exp(y);
      dy = ds_eff/sqrt(3+8*Q2);
      if (flag) {
        grid.templates[count].qual = Q;
        grid.templates[count].freq = f;
      }
      count++;
    }
  }
  return count;
}

void qn_template_grid(float dl, struct qnScope *grid)
{
  (*grid).n_tmplt = qn_set_grid(dl,*grid,0);
  if ((*grid).n_tmplt<1) error("qn_template_grid(): there are NO templates in this grid.");
  (*grid).templates =
    (struct qnTemplate*)malloc((size_t)((*grid).n_tmplt*sizeof(struct qnTemplate)));
  if (!((*grid).templates)) error("qn_template_grid(): could not allocate %d bytes memory.",
						((*grid).n_tmplt*sizeof(struct qnTemplate)));
  qn_set_grid(dl,*grid,1);
  return;
}

/* other routines */

void error(char *fmt, ...)
{
  va_list args;

  va_start(args, fmt);
  fprintf(stderr, "GRASP: error: ");
  vfprintf(stderr, fmt, args);
  fprintf(stderr,"\n%s\n",rcsid);
  va_end(args);
  exit(1);
}

int warning(char *fmt, ...)
{
  va_list args;

  va_start(args, fmt);
  fprintf(stderr, "GRASP: warning: ");
  vfprintf(stderr, fmt, args);
  fprintf(stderr,"\n%s\n",rcsid);
  va_end(args);
  return 1;
}
