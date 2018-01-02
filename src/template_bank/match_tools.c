/* GRASP: Copyright 1997,1998  Bruce Allen */
/* GRASP: Copyright 1997, Bruce Allen */
/* This routine: June/July 1997 Scott Hughes */

static char *rcsid="$Id: match_tools.c,v 1.2 1998/01/23 17:59:51 ballen Exp $\
\n$Name: RELEASE_1_9_8 $";

#include "grasp.h"
#define MAXFITPTS 300   /* Maximum number of points we keep for the fit */
#define HSCALE 1.e15    /* A convenient scaling factor; results independent
			   of it */
#define NCOF_cubic 7
#define NCOF_parab 3

/* These are prototypes of Numerical Recipes routines */
float **matrix(long,long,long,long);
void free_matrix(float **,long,long,long,long);
void gaussj(float **,int,float **,int);

int match_cubic(float m1ch,float m2ch,float matchcont,int order,
		float srate,float flo,float ftau,char *noise_file,
		float *semiaxis0,float *semiaxis1,float *theta,float mc[])
{
  void realft(float *,unsigned long,int);
  double mtot,eta,tau0save,tau1save,*fitpower;
  float *ch0tilde,*ch90tilde,freq,delta_f,*twice_inv_noise;
  float zcont,s_n0,s_n90,inverse_distance_scale,tc;
  float gammq(float,float),sh(float);
  float m1,m2,dtau0,dtau1;
  /* In the notation I've been using,
     x = tau0save - tau0,
     y = tau1save - tau1,
     z = match - 1 */
  float *x,*y,*z;
  /* things used in construction of ellipses */
  float *match_coef,aa,bb,cc,dd,ee,ff,gg,tmp;
  float eigval0,eigval1,*semiaxis,*eigvec0,*eigvec1,norm0,norm1;
  float phi,A0,A1,fit,fiterror=20;
  int fitindex=0,ellpts=25,ellpts_used;
  int i,j,k,okflag=0,npoint=131072,chirppoints=131072;

  /* Sanity check: if matchcont >=1, the functions's inputs
     are no good. */

  if(matchcont >=1 ) {
    GR_start_error("match_cubic()",rcsid,__FILE__,__LINE__);
    GR_report_error("The match contour value must be less than 1.\n");
    GR_end_error();
    exit(0);
  }
  zcont=matchcont-1.;

  x=(float *)malloc(sizeof(float)*(MAXFITPTS+ellpts+80));
  y=(float *)malloc(sizeof(float)*(MAXFITPTS+ellpts+80));
  z=(float *)malloc(sizeof(float)*(MAXFITPTS+ellpts+80));
  eigvec0=(float *)malloc(sizeof(float)*2);
  eigvec1=(float *)malloc(sizeof(float)*2);
  semiaxis=(float *)malloc(sizeof(float)*2);
  match_coef=(float *)malloc(sizeof(float)*7);
  
  /* Phase 0 and phase pi/2 chirps, in frequency domain.

     Keep increasing npoint by a factor of 2 until the number
     of points allocated is sufficient.

     POTENTIAL DANGER: If your computer doesn't have much
     RAM, this block of code is guaranteed to crash if you
     are trying to make chirps for small masses! */
  while (!okflag) {
    ch0tilde=(float *)malloc(sizeof(float)*npoint);
    ch90tilde=(float *)malloc(sizeof(float)*npoint);
    make_filters(m1ch,m2ch,ch0tilde,ch90tilde,flo,npoint,srate,
		 &chirppoints,&tc,0,order);
    if(chirppoints>npoint/2) {
      npoint*=2;
      free((void*)ch0tilde);
      free((void*)ch90tilde);
    } else okflag=1;
  }
  
  inverse_distance_scale=2.*HSCALE*(TSOLAR*C_LIGHT/MPC);
  for (i=0;i<chirppoints;i++){
    ch0tilde[i]*=inverse_distance_scale;
    ch90tilde[i]*=inverse_distance_scale;
  }  
  
  /* Twice the inverse of the mean noise power spectrum */
  twice_inv_noise=(float *)malloc(sizeof(float)*(npoint/2+1));
  fitpower=(double *)malloc(sizeof(double)*(npoint/2+1));
  
  /* and FFT the chirps */
  realft(ch0tilde-1,npoint,1);
  realft(ch90tilde-1,npoint,1);
  
  /* power spectrum */
  delta_f=srate/((float)npoint);
  noise_power(noise_file,npoint/2+1,delta_f,fitpower);
  for (i=0;i<npoint/2+1;i++) {
    freq = i*delta_f;
    fitpower[i]*=HSCALE*HSCALE;
    if (fitpower[i]!=0. && freq > flo)
      twice_inv_noise[i]=2./fitpower[i];
    else
      twice_inv_noise[i]=0.;
  }
  /* orthogonalize the chirps: we never modify ch0tilde, only ch90tilde */
  orthonormalize(ch0tilde,ch90tilde,twice_inv_noise,npoint,&s_n0,&s_n90);

  tau_of_mass(m1ch,m2ch,M_PI*ftau,&tau0save,&tau1save);

  /* For small values of the z contour (ie, as the match contour
     approaches 1), make the initial stepsize smaller. */
  dtau0 = (2.e-4)*sqrt(-zcont/.03);
  dtau1 = (5.e-5)*sqrt(-zcont/.03);

  /* Put together first grid of points, find initial guess for the
     equimatch contour ellipse.  This will give a lousy approximation
     to the ellipse, that I will then iterate on with a root finder. */
  k=0;
  while(fitindex < 9) {
    i = (k%3)-1; j = (k/3)-1;
    x[fitindex]= i*dtau0; y[fitindex]= j*dtau1;
    if(m_and_eta(tau0save+x[fitindex],tau1save+y[fitindex],&mtot,
		 &eta,1.e-6,1.e6,M_PI*ftau)) {
      m1 = 0.5*mtot*(1.-sqrt(1.-4.*eta)); m2 = mtot-m1;
      z[fitindex]=compute_match(m1,m2,ch0tilde,ch90tilde,
				inverse_distance_scale,twice_inv_noise,
				flo,s_n0,s_n90,npoint,srate,4000,order)-1.;
      /* Want to make sure that none of the points that are sampled
	 lie outside the eventual ellipse that will be found - such
	 data points can poison the entire sample.  */
      if(1.+z[fitindex] < .995*matchcont) {
	/* If the match found is much less than matchcont, the data
	   set is trash.  Cut dtau0,dtau1 down and start over. */
	dtau0 /= 3.;
	dtau1 /= 3.;
	k = -1;
	fitindex = -1;
      }
      fitindex++;
    }
    k++;
  }
  cubic(x-1,y-1,z-1,fitindex,match_coef-1);
  aa=match_coef[0]; bb=0.5*match_coef[1]; cc=match_coef[2];
  dd=match_coef[3]; ee=match_coef[4]; ff=match_coef[5]; gg=match_coef[6];
  /* Find eigenvalues and eigenvectors of the matrix {{aa,bb},{bb,cc}}.
     Using analytic formulae. */
  tmp = 0.5*(aa-cc)/bb;
  /* If the initial set of points is really lousy, the first "ellipse"
     won't even be elliptical: At least one of the eigenvalues will
     have the wrong sign.  Don't lose any sleep over it; we only need
     the eigenvalues to help decide how far we're going to step around
     when we try to find more points.  Just artificially change the
     eigenvalue's sign. */
  eigval0 = 0.5*(aa+cc) + fabs(bb)*sqrt(tmp*tmp+1.);
  if(eigval0 > 0.0) eigval0 *= -1.;
  eigval1 = 0.5*(aa+cc) - fabs(bb)*sqrt(tmp*tmp+1.);
  if(eigval1 > 0.0) eigval1 *= -1.;
  norm0 = sqrt(2.*(tmp*tmp + 1. + tmp*sqrt(tmp*tmp+1.)));
  norm1 = sqrt(2.*(tmp*tmp + 1. - tmp*sqrt(tmp*tmp+1.)));
  /* An important - and confusing - point: |eigval0| < |eigval1|,
     because the matrix is NEGATIVE DEFINITE.
     Therefore, semiaxis[0] > semiaxis[1]. */
  eigvec0[0] = (tmp + sqrt(tmp*tmp+1.))/norm0;
  eigvec1[0] = (tmp - sqrt(tmp*tmp+1.))/norm1;
  eigvec0[1] = 1./norm0;
  eigvec1[1] = 1./norm1;
  semiaxis[0] = sqrt(zcont/eigval0); semiaxis[1] = sqrt(zcont/eigval1);
  /* fiterror is a chi-squared-like statistic.  The fit is "good"
     if fiterror is in the neighborhood of unity.  This loop keeps
     gathering more data until the fit evaluated on the ellipse
     has fiterror less than 5. */
  while (fiterror > 5. && fitindex < (MAXFITPTS - 50)) {
    bisect(x,y,z,matchcont,&fitindex,eigvec0,eigvec1,srate,flo,ftau,
	   s_n0,s_n90,twice_inv_noise,inverse_distance_scale,
	   tau0save,tau1save,semiaxis,0,ch0tilde,ch90tilde,npoint);
    bisect(x,y,z,matchcont,&fitindex,eigvec0,eigvec1,srate,flo,ftau,
	   s_n0,s_n90,twice_inv_noise,inverse_distance_scale,
	   tau0save,tau1save,semiaxis,1,ch0tilde,ch90tilde,npoint);
    bisect(x,y,z,matchcont,&fitindex,eigvec0,eigvec1,srate,flo,ftau,
	   s_n0,s_n90,twice_inv_noise,inverse_distance_scale,
	   tau0save,tau1save,semiaxis,2,ch0tilde,ch90tilde,npoint);
    bisect(x,y,z,matchcont,&fitindex,eigvec0,eigvec1,srate,flo,ftau,
	   s_n0,s_n90,twice_inv_noise,inverse_distance_scale,
	   tau0save,tau1save,semiaxis,3,ch0tilde,ch90tilde,npoint);
    cubic(x-1,y-1,z-1,fitindex,match_coef-1);
    aa=match_coef[0]; bb=0.5*match_coef[1]; cc=match_coef[2];
    dd=match_coef[3]; ee=match_coef[4]; ff=match_coef[5]; gg=match_coef[6];
    tmp = 0.5*(aa-cc)/bb;
    eigval0 = 0.5*(aa+cc) + fabs(bb)*sqrt(tmp*tmp+1.);
    if(eigval0 > 0.0) eigval0 *= -1.;
    eigval1 = 0.5*(aa+cc) - fabs(bb)*sqrt(tmp*tmp+1.);
    if(eigval1 > 0.0) eigval1 *= -1.;
    norm0 = sqrt(2.*(tmp*tmp + 1. + tmp*sqrt(tmp*tmp+1.)));
    norm1 = sqrt(2.*(tmp*tmp + 1. - tmp*sqrt(tmp*tmp+1.)));
    eigvec0[0] = (tmp + sqrt(tmp*tmp+1.))/norm0;
    eigvec1[0] = (tmp - sqrt(tmp*tmp+1.))/norm1;
    eigvec0[1] = 1./norm0;
    eigvec1[1] = 1./norm1;
    semiaxis[0] = sqrt(zcont/eigval0); semiaxis[1] = sqrt(zcont/eigval1);
    /* Now, step around the ellipse that has been found.  Code will
       use at most ellpts points (possibly less if a piece of the
       ellipse lies south of the equal mass line).  Accumulate the
       fiterror on this ellipse. */
    fiterror=0.0; ellpts_used=0;
    for (i=0;i<ellpts;i++) {
      phi = i*2.0*M_PI/((float)ellpts);
      A0 = semiaxis[0]*cos(phi); A1 = semiaxis[1]*sin(phi);
      x[fitindex] = A0*eigvec0[0] + A1*eigvec1[0];
      y[fitindex] = A0*eigvec0[1] + A1*eigvec1[1];
      if(m_and_eta(tau0save+x[fitindex],tau1save+y[fitindex],
		   &mtot,&eta,1.e-6,1.e6,M_PI*ftau)) {
	m1 = 0.5*mtot*(1.-sqrt(1.-4.*eta)); m2 = mtot-m1;
	z[fitindex]=compute_match(m1,m2,ch0tilde,ch90tilde,
				  inverse_distance_scale,twice_inv_noise,
				  flo,s_n0,s_n90,npoint,srate,4000,4)-1.;
	fit = 1.+aa*x[fitindex]*x[fitindex]+2.0*bb*x[fitindex]*y[fitindex]
	  +cc*y[fitindex]*y[fitindex]+dd*x[fitindex]*x[fitindex]*x[fitindex]
	  +ee*y[fitindex]*y[fitindex]*y[fitindex]
	  +ff*x[fitindex]*x[fitindex]*y[fitindex]
	  +gg*x[fitindex]*y[fitindex]*y[fitindex];
	fiterror += (fit-(z[fitindex]+1.))*(fit-(z[fitindex]+1.))*1.e6;
	fitindex++; ellpts_used++;
      }
    }
    fiterror /= ((float)ellpts_used);
  }
  *semiaxis0=semiaxis[0]; *semiaxis1=semiaxis[1];
  *theta=atan2(eigvec0[1],eigvec0[0]);
  mc[0]=aa; mc[1]=2.0*bb; mc[2]=cc; mc[3]=dd; mc[4]=ee;
  mc[5]=ff; mc[6]=gg;
  if (fitindex >= (MAXFITPTS - 50)) {
    /* Code did not find a decent fit. */
    free((void*)(x));
    free((void*)(y));
    free((void*)(z));
    free((void*)(eigvec0));
    free((void*)(eigvec1));
    free((void*)(semiaxis));
    free((void*)(match_coef));
    free((void*)(ch0tilde));
    free((void*)(ch90tilde));
    free((void*)(twice_inv_noise));
    free((void*)(fitpower));
    return(0);
  } else {
    free((void*)(x));
    free((void*)(y));
    free((void*)(z));
    free((void*)(eigvec0));
    free((void*)(eigvec1));
    free((void*)(semiaxis));
    free((void*)(match_coef));
    free((void*)(ch0tilde));
    free((void*)(ch90tilde));
    free((void*)(twice_inv_noise));
    free((void*)(fitpower));
    return(1);
  }
}

int match_parab(float m1ch,float m2ch,float matchcont,int order,
		float srate,float flo,float ftau,char *noise_file,
		float *semiaxis0,float *semiaxis1,float *theta,float mc[])
{
  void realft(float *,unsigned long,int);
  double mtot,eta,tau0save,tau1save,*fitpower;
  float *ch0tilde,*ch90tilde,freq,delta_f,*twice_inv_noise;
  float zcont,s_n0,s_n90,inverse_distance_scale,tc;
  float gammq(float,float),sh(float);
  float m1,m2,dtau0,dtau1;
  /* In the notation I've been using,
     x = tau0save - tau0,
     y = tau1save - tau1,
     z = match - 1 */
  float *x,*y,*z;
  /* coefficients used in polynomial fits */
  float *match_coef,aa,bb,cc,tmp;
  /* things used in construction of ellipses */
  float eigval0,eigval1,*semiaxis,*eigvec0,*eigvec1,norm0,norm1;
  float phi,A0,A1,fit,fiterror=20;
  int fitindex=0,ellpts=25,ellpts_used;
  int i,j,k,okflag=0,npoint=131072,chirppoints=131072;

  /* Sanity check: if matchcont >=1, the functions's inputs
     are no good. */

  if(matchcont >=1 ) {
    GR_start_error("match_cubic()",rcsid,__FILE__,__LINE__);
    GR_report_error("The match contour value must be less than 1.\n");
    GR_end_error();
    exit(0);
  }
  zcont=matchcont-1.;

  x=(float *)malloc(sizeof(float)*(MAXFITPTS+ellpts+80));
  y=(float *)malloc(sizeof(float)*(MAXFITPTS+ellpts+80));
  z=(float *)malloc(sizeof(float)*(MAXFITPTS+ellpts+80));
  eigvec0=(float *)malloc(sizeof(float)*2);
  eigvec1=(float *)malloc(sizeof(float)*2);
  semiaxis=(float *)malloc(sizeof(float)*2);
  match_coef=(float *)malloc(sizeof(float)*3);
  
  /* Phase 0 and phase pi/2 chirps, in frequency domain.

     Keep increasing npoint by a factor of 2 until the number
     of points allocated is sufficient.

     POTENTIAL DANGER: If your computer doesn't have much
     RAM, this block of code is guaranteed to crash if you
     are trying to make chirps for small masses! */
  while (!okflag) {
    ch0tilde=(float *)malloc(sizeof(float)*npoint);
    ch90tilde=(float *)malloc(sizeof(float)*npoint);
    make_filters(m1ch,m2ch,ch0tilde,ch90tilde,flo,npoint,srate,
		 &chirppoints,&tc,0,order);
    if(chirppoints>npoint/2) {
      npoint*=2;
      free((void*)ch0tilde);
      free((void*)ch90tilde);
    } else okflag=1;
  }
  
  inverse_distance_scale=2.*HSCALE*(TSOLAR*C_LIGHT/MPC);
  for (i=0;i<chirppoints;i++){
    ch0tilde[i]*=inverse_distance_scale;
    ch90tilde[i]*=inverse_distance_scale;
  }  
  
  /* Twice the inverse of the mean noise power spectrum */
  twice_inv_noise=(float *)malloc(sizeof(float)*(npoint/2+1));
  fitpower=(double *)malloc(sizeof(double)*(npoint/2+1));
  
  /* and FFT the chirps */
  realft(ch0tilde-1,npoint,1);
  realft(ch90tilde-1,npoint,1);

  /* power spectrum */
  delta_f=srate/((float)npoint);
  noise_power(noise_file,npoint/2+1,delta_f,fitpower);
  for (i=0;i<npoint/2+1;i++) {
    freq = i*delta_f;
    fitpower[i]*=HSCALE*HSCALE;
    if (fitpower[i]!=0. && freq > flo)
      twice_inv_noise[i]=2./fitpower[i];
    else
      twice_inv_noise[i]=0.;
  }

  /* orthogonalize the chirps: we never modify ch0tilde, only ch90tilde */
  orthonormalize(ch0tilde,ch90tilde,twice_inv_noise,npoint,&s_n0,&s_n90);

  tau_of_mass(m1ch,m2ch,M_PI*ftau,&tau0save,&tau1save);

  /* For small values of the z contour (ie, as the match contour
     approaches 1), make the initial stepsize smaller. */
  dtau0 = (2.e-4)*sqrt(-zcont/.03);
  dtau1 = (5.e-5)*sqrt(-zcont/.03);

  /* Put together first grid of points, find initial guess for the
     equimatch contour ellipse.  This will give a lousy approximation
     to the ellipse, that I will then iterate on with a root finder. */
  k=0;
  while(fitindex < 9) {
    i = (k%3)-1; j = (k/3)-1;
    x[fitindex]= i*dtau0; y[fitindex]= j*dtau1;
    if(m_and_eta(tau0save+x[fitindex],tau1save+y[fitindex],&mtot,
		 &eta,1.e-6,1.e6,M_PI*ftau)) {
      m1 = 0.5*mtot*(1.-sqrt(1.-4.*eta)); m2 = mtot-m1;
      z[fitindex]=compute_match(m1,m2,ch0tilde,ch90tilde,
				inverse_distance_scale,twice_inv_noise,
				flo,s_n0,s_n90,npoint,srate,4000,order)-1.;
      /* Want to make sure that none of the points that are sampled
	 lie outside the eventual ellipse that will be found - such
	 data points can poison the entire sample.  */
      if(1.+z[fitindex] < .995*matchcont) {
	/* If the match found is much less than matchcont, the data
	   set is trash.  Cut dtau0,dtau1 down and start over. */
	dtau0 /= 3.;
	dtau1 /= 3.;
	k = -1;
	fitindex = -1;
      }
      fitindex++;
    }
    k++;
  }
  parab(x-1,y-1,z-1,fitindex,match_coef-1);
  aa=match_coef[0]; bb=0.5*match_coef[1]; cc=match_coef[2];
  /* Find eigenvalues and eigenvectors of the matrix {{aa,bb},{bb,cc}}.
     Using analytic formulae. */
  tmp = 0.5*(aa-cc)/bb;
  /* If the initial set of points is really lousy, the first "ellipse"
     won't even be elliptical: At least one of the eigenvalues will
     have the wrong sign.  Don't lose any sleep over it; we only need
     the eigenvalues to help decide how far we're going to step around
     when we try to find more points.  Just artificially change the
     eigenvalue's sign. */
  eigval0 = 0.5*(aa+cc) + fabs(bb)*sqrt(tmp*tmp+1.);
  if(eigval0 > 0.0) eigval0 *= -1.;
  eigval1 = 0.5*(aa+cc) - fabs(bb)*sqrt(tmp*tmp+1.);
  if(eigval1 > 0.0) eigval1 *= -1.;
  norm0 = sqrt(2.*(tmp*tmp + 1. + tmp*sqrt(tmp*tmp+1.)));
  norm1 = sqrt(2.*(tmp*tmp + 1. - tmp*sqrt(tmp*tmp+1.)));
  /* An important - and confusing - point: |eigval0| < |eigval1|,
     because the matrix is NEGATIVE DEFINITE.
     Therefore, semiaxis[0] > semiaxis[1]. */
  eigvec0[0] = (tmp + sqrt(tmp*tmp+1.))/norm0;
  eigvec1[0] = (tmp - sqrt(tmp*tmp+1.))/norm1;
  eigvec0[1] = 1./norm0;
  eigvec1[1] = 1./norm1;
  semiaxis[0] = sqrt(zcont/eigval0); semiaxis[1] = sqrt(zcont/eigval1);
  /* fiterror is a chi-squared-like statistic.  The fit is "good"
     if fiterror is in the neighborhood of unity.  This loop keeps
     gathering more data until the fit evaluated on the ellipse
     has fiterror less than 5. */
  while (fiterror > 5. && fitindex < (MAXFITPTS - 50)) {
    bisect(x,y,z,matchcont,&fitindex,eigvec0,eigvec1,srate,flo,ftau,
	   s_n0,s_n90,twice_inv_noise,inverse_distance_scale,
	   tau0save,tau1save,semiaxis,0,ch0tilde,ch90tilde,npoint);
    bisect(x,y,z,matchcont,&fitindex,eigvec0,eigvec1,srate,flo,ftau,
	   s_n0,s_n90,twice_inv_noise,inverse_distance_scale,
	   tau0save,tau1save,semiaxis,1,ch0tilde,ch90tilde,npoint);
    bisect(x,y,z,matchcont,&fitindex,eigvec0,eigvec1,srate,flo,ftau,
	   s_n0,s_n90,twice_inv_noise,inverse_distance_scale,
	   tau0save,tau1save,semiaxis,2,ch0tilde,ch90tilde,npoint);
    bisect(x,y,z,matchcont,&fitindex,eigvec0,eigvec1,srate,flo,ftau,
	   s_n0,s_n90,twice_inv_noise,inverse_distance_scale,
	   tau0save,tau1save,semiaxis,3,ch0tilde,ch90tilde,npoint);
    parab(x-1,y-1,z-1,fitindex,match_coef-1);
    aa=match_coef[0]; bb=0.5*match_coef[1]; cc=match_coef[2];
    tmp = 0.5*(aa-cc)/bb;
    eigval0 = 0.5*(aa+cc) + fabs(bb)*sqrt(tmp*tmp+1.);
    if(eigval0 > 0.0) eigval0 *= -1.;
    eigval1 = 0.5*(aa+cc) - fabs(bb)*sqrt(tmp*tmp+1.);
    if(eigval1 > 0.0) eigval1 *= -1.;
    norm0 = sqrt(2.*(tmp*tmp + 1. + tmp*sqrt(tmp*tmp+1.)));
    norm1 = sqrt(2.*(tmp*tmp + 1. - tmp*sqrt(tmp*tmp+1.)));
    eigvec0[0] = (tmp + sqrt(tmp*tmp+1.))/norm0;
    eigvec1[0] = (tmp - sqrt(tmp*tmp+1.))/norm1;
    eigvec0[1] = 1./norm0;
    eigvec1[1] = 1./norm1;
    semiaxis[0] = sqrt(zcont/eigval0); semiaxis[1] = sqrt(zcont/eigval1);
    /* Now, step around the ellipse that has been found.  Code will
       use at most ellpts points (possibly less if a piece of the
       ellipse lies south of the equal mass line).  Accumulate the
       fiterror on this ellipse. */
    fiterror=0.0; ellpts_used=0;
    for (i=0;i<ellpts;i++) {
      phi = i*2.0*M_PI/((float)ellpts);
      A0 = semiaxis[0]*cos(phi); A1 = semiaxis[1]*sin(phi);
      x[fitindex] = A0*eigvec0[0] + A1*eigvec1[0];
      y[fitindex] = A0*eigvec0[1] + A1*eigvec1[1];
      if(m_and_eta(tau0save+x[fitindex],tau1save+y[fitindex],
		   &mtot,&eta,1.e-6,1.e6,M_PI*ftau)) {
	m1 = 0.5*mtot*(1.-sqrt(1.-4.*eta)); m2 = mtot-m1;
	z[fitindex]=compute_match(m1,m2,ch0tilde,ch90tilde,
				  inverse_distance_scale,twice_inv_noise,flo,
				  s_n0,s_n90,npoint,srate,4000,4)-1.;
	fit = 1.+aa*x[fitindex]*x[fitindex]+2.0*bb*x[fitindex]*y[fitindex]
	  +cc*y[fitindex]*y[fitindex];
	fiterror += (fit-(z[fitindex]+1.))*(fit-(z[fitindex]+1.))*1.e6;
	fitindex++; ellpts_used++;
      }
    }
    fiterror /= ((float)ellpts_used);
  }
  *semiaxis0=semiaxis[0]; *semiaxis1=semiaxis[1];
  *theta=atan2(eigvec0[1],eigvec0[0]);
  mc[0]=aa; mc[1]=2.0*bb; mc[2]=cc;
  if (fitindex >= (MAXFITPTS-50)) {
    /* Code did not find a decent fit. */
    free((void*)(x));
    free((void*)(y));
    free((void*)(z));
    free((void*)(eigvec0));
    free((void*)(eigvec1));
    free((void*)(semiaxis));
    free((void*)(match_coef));
    free((void*)(ch0tilde));
    free((void*)(ch90tilde));
    free((void*)(twice_inv_noise));
    free((void*)(fitpower));
    return(0);
  } else {
    free((void*)(x));
    free((void*)(y));
    free((void*)(z));
    free((void*)(eigvec0));
    free((void*)(eigvec1));
    free((void*)(semiaxis));
    free((void*)(match_coef));
    free((void*)(ch0tilde));
    free((void*)(ch90tilde));
    free((void*)(twice_inv_noise));
    free((void*)(fitpower));
    return(1);
  }
}

void bisect(float x[],float y[],float z[],float matchcont,int *index,
	    float eigvec0[],float eigvec1[],float srate,float flo,
	    float ftau,float s_n0,float s_n90,float twice_inv_noise[],
	    float inverse_distance_scale,float tau0,float tau1,
	    float semiaxis[],int bisectflag,float ch0tilde[],
	    float ch90tilde[],int npoint)
{
  void realft(float *,unsigned long,int);
  float *A,*dA,*oldA;
  float stepsign;
  double mtot,eta;
  float m1,m2;
  float match=1;
  int stepflag,holdflag;

  /* If bisectflag==0, bisect along positive A0 direction.
     If bisectflag==1, bisect along negative A0 direction.
     If bisectflag==2, bisect along positive A1 direction.
     If bisectflag==3, bisect along negative A1 direction. */
  if (bisectflag < 0 || bisectflag > 3) {
    GR_start_error("bisect()",rcsid,__FILE__,__LINE__);
    GR_report_error("Invalid bisectflag %d in function bisect\n",bisectflag);
    GR_end_error();
    exit(0);
  }
  A=(float *)malloc(sizeof(float)*2);
  dA=(float *)malloc(sizeof(float)*2);
  oldA=(float *)malloc(sizeof(float)*2);
  
  if(bisectflag<2) {
    stepflag=0; holdflag=1;
  } else {
    stepflag=1; holdflag=0;
  }
  if(bisectflag%2) stepsign= -1.;
  else stepsign=1.;

  A[stepflag] = oldA[stepflag] = stepsign*semiaxis[stepflag];
  A[holdflag] = oldA[holdflag] = 0.0;
  dA[stepflag] = .4*A[stepflag];
  dA[holdflag] = 0.;

  /* First, check to see if the end of the axis is even an allowed
     point.  If it isn't, forget it: the routine will get trapped. */
  x[*index] = A[0]*eigvec0[0] + A[1]*eigvec1[0];
  y[*index] = A[0]*eigvec0[1] + A[1]*eigvec1[1];
  if(!m_and_eta(tau0+x[*index],tau1+y[*index],&mtot,&eta,
	       1.e-6,1.e6,M_PI*ftau))
    return;

  while(fabs(match-matchcont) > .0002 &&
	fabs(dA[stepflag]/A[stepflag]) > 5.e-4) {
    A[0] += dA[0]; A[1] += dA[1];
    x[*index] = A[0]*eigvec0[0] + A[1]*eigvec1[0];
    y[*index] = A[0]*eigvec0[1] + A[1]*eigvec1[1];
    if(m_and_eta(tau0+x[*index],tau1+y[*index],&mtot,&eta,
		 1.e-6,1.e6,M_PI*ftau)) {
      m1 = 0.5*mtot*(1.-sqrt(1.-4.*eta)); m2 = mtot-m1;
      match=compute_match(m1,m2,ch0tilde,ch90tilde,inverse_distance_scale,
			  twice_inv_noise,flo,s_n0,s_n90,npoint,srate,
			  4000,4);
      z[*index]=match-1.0;
      if (match > matchcont) {
	if (fabs(oldA[stepflag]) > fabs(A[stepflag])) {
	  /* Previous change of dA[stepflag] occurred on the other
	     side of the root: time to change.  Note strange
	     reduction factor; this is done so that it doesn't (often)
	     step on the same point twice. */
	  dA[stepflag] *= -M_PI/10.0;
	  oldA[stepflag] = A[stepflag];
	}
      } else {
	if (fabs(oldA[stepflag]) < fabs(A[stepflag])) {
	  /* Previous change of dA[stepflag] occurred on the other
	     side of the root: time to change.  Note strange
	     reduction factor; this is done so that it doesn't (often)
	     step on the same point twice. */
	  dA[stepflag] *= -M_PI/10.0;
	  oldA[stepflag] = A[stepflag];
	}
      }
      (*index)++;
    } else {
      /* Stepped into the land of serpents and dragons, where masses
	 exist not. Back out and reduce the stepsize. */
      A[stepflag] -= dA[stepflag];
      dA[stepflag] *= 0.5;
    }
  }
  free((void*)(A));
  free((void*)(dA));
  free((void*)(oldA));
  return;
}

/* Note: assumes arrays are UNIT offset */

int cubic(float x[], float y[], float z[], int Ndat, float a[])
{
  int i,j,k;
  float **beta,**alpha,
    Xfunc_cubic(float x[],float y[],int i,int k);

  beta = matrix(1,NCOF_cubic,1,1);
  alpha = matrix(1,NCOF_cubic,1,NCOF_cubic);
  for(j=1;j<=NCOF_cubic;j++) {
    for(k=1;k<=NCOF_cubic;k++) {
      for(i=1;i<=Ndat;i++) {
	if(i==1) {
	  if(j==1) beta[k][1] = 0.0;
	  alpha[k][j] = 0.0;
	}
	if(j==1)
	  beta[k][1] += z[i]*Xfunc_cubic(x,y,i,k);
	alpha[k][j] += Xfunc_cubic(x,y,i,k)*Xfunc_cubic(x,y,i,j);
      }
    }
  }
  gaussj(alpha,NCOF_cubic,beta,1);
  for(i=1;i<=NCOF_cubic;i++) a[i]=beta[i][1];
  free_matrix(alpha,1,NCOF_cubic,1,NCOF_cubic);
  free_matrix(beta,1,NCOF_cubic,1,1);
  return 1;
}

/* In the notation of my notes, returns $X_k(x_i,y_i)$ */
float Xfunc_cubic(float x[], float y[], int i, int k)
{
  if(k==1)
    return x[i]*x[i];
  else if(k==2)
    return x[i]*y[i];
  else if(k==3)
    return y[i]*y[i];
  else if(k==4)
    return x[i]*x[i]*x[i];
  else if(k==5)
    return y[i]*y[i]*y[i];
  else if(k==6)
    return x[i]*x[i]*y[i];
  else return x[i]*y[i]*y[i];
}

int parab(float x[], float y[], float z[], int Ndat, float a[])
{
  int i,j,k;
  float **beta,**alpha,
    Xfunc_parab(float x[],float y[],int i,int k);

  beta = matrix(1,NCOF_parab,1,1);
  alpha = matrix(1,NCOF_parab,1,NCOF_parab);
  for(j=1;j<=NCOF_parab;j++) {
    for(k=1;k<=NCOF_parab;k++) {
      for(i=1;i<=Ndat;i++) {
	if(i==1) {
	  if(j==1) beta[k][1] = 0.0;
	  alpha[k][j] = 0.0;
	}
	if(j==1)
	  beta[k][1] += z[i]*Xfunc_parab(x,y,i,k);
	alpha[k][j] += Xfunc_parab(x,y,i,k)*Xfunc_parab(x,y,i,j);
      }
    }
  }
  gaussj(alpha,NCOF_parab,beta,1);
  for(i=1;i<=NCOF_parab;i++) a[i]=beta[i][1];
  free_matrix(alpha,1,NCOF_parab,1,NCOF_parab);
  free_matrix(beta,1,NCOF_parab,1,1);
  return 1;
}

/* In the notation of my notes, returns $X_k(x_i,y_i)$ */
float Xfunc_parab(float x[], float y[], int i, int k)
{
  if(k==1)
    return x[i]*x[i];
  else if(k==2)
    return x[i]*y[i];
  else
    return y[i]*y[i];
}

float compute_match(float m1,float m2,float ch0tilde[],float ch90tilde[],
		float inverse_distance_scale,float twice_inv_noise[],
		float flo,float s_n0,float s_n90,int npoint,float srate,
		int err_cd_sprs,int order)
{
  void realft(float *,unsigned long,int);
  int i,chirppoints,maxi;
  float *output0,*output90,*template0tilde,*template90tilde;
  float tc,t_n0,t_n90,lin0,lin90,snr_max,var;

  output0=(float *)malloc(sizeof(float)*npoint);
  output90=(float *)malloc(sizeof(float)*npoint);
  template0tilde=(float *)malloc(sizeof(float)*npoint);
  template90tilde=(float *)malloc(sizeof(float)*npoint);
  
  make_filters(m1,m2,template0tilde,template90tilde,flo,npoint,
	       srate,&chirppoints,&tc,err_cd_sprs,order);
  for(i=0;i<npoint;i++) {
    template0tilde[i]*=inverse_distance_scale;
    template90tilde[i]*=inverse_distance_scale;
  }
  realft(template0tilde-1,npoint,1);
  realft(template90tilde-1,npoint,1);
  orthonormalize(template0tilde,template90tilde,twice_inv_noise,
		 npoint,&t_n0,&t_n90);
  find_chirp(template0tilde,ch0tilde,ch90tilde,twice_inv_noise,
	     s_n0,s_n90,output0,output90,npoint,0,&maxi,
	     &snr_max,&lin0,&lin90,&var);
  free((void*)(template0tilde));
  free((void*)(template90tilde));
  free((void*)(output0));
  free((void*)(output90));
  return(t_n0*snr_max/sqrt(2.0));
}

