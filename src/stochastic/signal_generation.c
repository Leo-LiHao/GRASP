/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"
static char *rcsid="$Id: signal_generation.c,v 1.6 1998/01/23 17:59:49 ballen Exp $\
\n$Name: RELEASE_1_9_8 $";

float gasdev(int *);
void four1(float *,int,int);
void realft(float *,int,int);

/* Function name: monte_carlo() */

void monte_carlo(int fake_sb,int fake_noise1,int fake_noise2,int n,
		 float delta_t,float omega_0,float f_low,float f_high,
		 double *gamma12,double *power1,double *power2,
		 double *whiten1,double *whiten2,
		 float *out1,float *out2,int *pseed)
{
  int i;

  static int last_n=0;
  static float *left1,*left2,*right1,*right2;
  static float *left_sb1,*left_sb2,*right_sb1,*right_sb2;
  static float *left_noise1,*left_noise2,*right_noise1,*right_noise2;

  if (n!=last_n) {
    
    /* (re)allocate memory if current n differs from previous n */
    left1=(float *)realloc(left1,n*sizeof(float));
    if (left1==NULL) {
      GR_start_error("monte_carlo()",rcsid,__FILE__,__LINE__);
      GR_report_error("Failed to allocate %d floats.\n",n);
	GR_end_error();
      abort();
    }
    left2=(float *)realloc(left2,n*sizeof(float));
    if (left2==NULL) {
      GR_start_error("monte_carlo()",rcsid,__FILE__,__LINE__);
      GR_report_error("Failed to allocate %d floats.\n",n);
	GR_end_error();
      abort();
    }
    right1=(float *)realloc(right1,n*sizeof(float));
    if (right1==NULL) {
      GR_start_error("monte_carlo()",rcsid,__FILE__,__LINE__);
      GR_report_error("Failed to allocate %d floats.\n",n);
	GR_end_error();
      abort();
    }
    right2=(float *)realloc(right2,n*sizeof(float));
    if (right2==NULL) {
      GR_start_error("monte_carlo()",rcsid,__FILE__,__LINE__);
      GR_report_error("Failed to allocate %d floats.\n",n);
	GR_end_error();
      abort();
    }
    left_sb1=(float *)realloc(left_sb1,n*sizeof(float));
    if (left_sb1==NULL) {
      GR_start_error("monte_carlo()",rcsid,__FILE__,__LINE__);
      GR_report_error("Failed to allocate %d floats.\n",n);
	GR_end_error();
      abort();
    }
    left_sb2=(float *)realloc(left_sb2,n*sizeof(float));
    if (left_sb2==NULL) {
      GR_start_error("monte_carlo()",rcsid,__FILE__,__LINE__);
      GR_report_error("Failed to allocate %d floats.\n",n);
	GR_end_error();
      abort();
    }
    right_sb1=(float *)realloc(right_sb1,n*sizeof(float));
    if (right_sb1==NULL) {
      GR_start_error("monte_carlo()",rcsid,__FILE__,__LINE__);
      GR_report_error("Failed to allocate %d floats.\n",n);
	GR_end_error();
      abort();
    }
    right_sb2=(float *)realloc(right_sb2,n*sizeof(float));
    if (right_sb2==NULL) {
      GR_start_error("monte_carlo()",rcsid,__FILE__,__LINE__);
      GR_report_error("Failed to allocate %d floats.\n",n);
	GR_end_error();
      abort();
    }
    left_noise1=(float *)realloc(left_noise1,n*sizeof(float));
    if (left_noise1==NULL) {
      GR_start_error("monte_carlo()",rcsid,__FILE__,__LINE__);
      GR_report_error("Failed to allocate %d floats.\n",n);
	GR_end_error();
      abort();
    }
    left_noise2=(float *)realloc(left_noise2,n*sizeof(float));
    if (left_noise2==NULL) {
      GR_start_error("monte_carlo()",rcsid,__FILE__,__LINE__);
      GR_report_error("Failed to allocate %d floats.\n",n);
	GR_end_error();
      abort();
    }
    right_noise1=(float *)realloc(right_noise1,n*sizeof(float));
    if (right_noise1==NULL) {
      GR_start_error("monte_carlo()",rcsid,__FILE__,__LINE__);
      GR_report_error("Failed to allocate %d floats.\n",n);
	GR_end_error();
      abort();
    }
    right_noise2=(float *)realloc(right_noise2,n*sizeof(float));
    if (right_noise2==NULL) {
      GR_start_error("monte_carlo()",rcsid,__FILE__,__LINE__);
      GR_report_error("Failed to allocate %d floats.\n",n);
	GR_end_error();
      abort();
    }
				  
    /* clear stochastic background and detector noise arrays */
    for (i=0;i<n;i++) {
      left_sb1[i]=left_sb2[i]=right_sb1[i]=right_sb2[i]=0.0;
      left_noise1[i]=left_noise2[i]=right_noise1[i]=right_noise2[i]=0.0;
    }
				  
    /* save value of n */
    last_n=n;

  } /* end if (n!=last_n) */

  /* simulate detector noise, if desired */
  if (fake_noise1==1) { /* site 1 */
    simulate_noise(n,delta_t,power1,whiten1,left_noise1,pseed);
    simulate_noise(n,delta_t,power1,whiten1,right_noise1,pseed);
  }
  if (fake_noise2==1) { /* site 2 */
    simulate_noise(n,delta_t,power2,whiten2,left_noise2,pseed);
    simulate_noise(n,delta_t,power2,whiten2,right_noise2,pseed);
  }

  /* simulate stochastic background, if desired */
  if (fake_sb==1) {
    simulate_sb(n,delta_t,omega_0,f_low,f_high,gamma12,whiten1,whiten2,
                left_sb1,left_sb2,pseed);
    simulate_sb(n,delta_t,omega_0,f_low,f_high,gamma12,whiten1,whiten2,
                right_sb1,right_sb2,pseed);
  }

  /* combine stochastic background and detector noise */
  for (i=0;i<n;i++) {
    left1[i]=left_sb1[i]+left_noise1[i];
    right1[i]=right_sb1[i]+right_noise1[i];
    left2[i]=left_sb2[i]+left_noise2[i];
    right2[i]=right_sb2[i]+right_noise2[i];
  }

  /* combine "left" and "right" data to get continuous-in-time output */
  combine_data(1,n,left1,right1,out1);
  combine_data(2,n,left2,right2,out2);

  return;
}
/* Filename: get_IFO12() */

void get_IFO12(FILE *fp1,FILE *fp2,FILE *fp1lock,FILE *fp2lock,
	       int n,float *out1,float *out2,float *srate1,float *srate2)
{
  int i,code1=0,code2=0,remain1,remain2;
  float tstart;

  static int last_n=0;
  static short *datas1,*datas2;

  /* (re)allocate memory if current n differs from previous n */

  if (n!=last_n) {

    datas1=(short *)realloc(datas1,n*sizeof(short));
    if (datas1==NULL) {
      GR_start_error("get_IFO12()",rcsid,__FILE__,__LINE__);
      GR_report_error("Failed to allocate %d shorts.\n",n);
	GR_end_error();
      abort();
    }
    datas2=(short *)realloc(datas2,n*sizeof(short));
    if (datas2==NULL) {
      GR_start_error("get_IFO12()",rcsid,__FILE__,__LINE__);
      GR_report_error("Failed to allocate %d shorts.\n",n);
	GR_end_error();
      abort();
    }

    /* save value of n */
    last_n=n;
  }

  /* get interferometer data */

  /* NOTE1:  we should skip 3 minutes into locked sections! */
  /* NOTE2:  if the data stream is interrupted, should we do something? */
 
  code1= get_data(fp1,fp1lock,&tstart,n,datas1,&remain1,srate1,0);
  code2=get_data2(fp2,fp2lock,&tstart,n,datas2,&remain2,srate2,0);

  /* put data into floats */
  for (i=0;i<n;i++) {
    out1[i]=(float)datas1[i];
    out2[i]=(float)datas2[i];
  }

  /* print warning message and clear output arrays if no more data */
  if (code1==0 || code2==0) {
    GR_start_error("monte_carlo()",rcsid,__FILE__,__LINE__);
    GR_report_error("No data remaining in one or both channels.\n");
    GR_report_error("No longer adding IFO data to output.\n");
	GR_end_error();

    /* set output arrays to zero */
    for (i=0;i<n;i++) out1[i]=out2[i]=0.0;

  }

  return;
}
/* Function name: simulate_sb() */

void simulate_sb(int n,float delta_t,float omega_0,float f_low,float f_high,
		 double *gamma12,double *whiten1,double *whiten2,
		 float *out1,float *out2,int *pseed)
{
  int i,pos,neg;
  float f,norm,prop,re1,im1,re2,im2;
  double re_w1,im_w1,re_w2,im_w2,real,imag;
  double g,g2,s;

  static int last_n=0;
  static float *data;

  /* (re)allocate memory if current n differs from previous n */

  if (n!=last_n) {
    data=(float *)realloc(data,2*n*sizeof(float));
    if (data==NULL) {
      GR_start_error("simulate_sb()",rcsid,__FILE__,__LINE__);
      GR_report_error("Failed to allocate %d floats.\n",2*n);
	GR_end_error();
      abort();
    }

    /* save value of n */
    last_n=n;
  }

  /* simulate stochastic background */

  /* frequency independent part of the normalization */
  norm=1.0/((float)n); /* for fft */
  norm*=sqrt((double)n/delta_t)*sqrt(3.0/20.0)*(HUBBLE/M_PI);
  norm*=sqrt((double)omega_0);

  /* positive frequencies */
  for (i=1;i<n/2;i++) {
    f=((float)i)/(n*delta_t);

    if (f<f_low || f>f_high) {
      im1=re1=im2=re2=0.0;
    }
    else {
      /* frequency dependence is f^-(3/2) */
      prop=norm/(f*sqrt((double)f));

      /* find real and imaginary parts, with random phases. */
      /* central limit theorem implies gaussian distribution */
      prop/=sqrt(2.0);
      re1=prop*gasdev(pseed);
      im1=prop*gasdev(pseed);

      /* fourier components at site 2 are correlated with those */
      /* at site 1 via the overlap reduction function           */
      g=gamma12[i];
      g2=g*g;
      s=sqrt(1.0-g2);
      re2=re1*g+s*prop*gasdev(pseed);
      im2=im1*g+s*prop*gasdev(pseed);
    }

    /* whiten data (site 1) */
    re_w1=whiten1[2*i];
    im_w1=whiten1[2*i+1];

    real=re1*re_w1-im1*im_w1;
    imag=re1*im_w1+im1*re_w1;
    re1=real;
    im1=imag;

    /* whiten data (site 2) */
    re_w2=whiten2[2*i];
    im_w2=whiten2[2*i+1];

    real=re2*re_w2-im2*im_w2;
    imag=re2*im_w2+im2*re_w2;
    re2=real;
    im2=imag;

    /* pack fourier components into data array */
    pos=2*i;
    neg=2*(n-i);

    /* real and imaginary parts, positive frequency */
    data[pos]=re1-im2;
    data[pos+1]=im1+re2;

    /* real and imaginary parts, negative frequency  */
    data[neg]=re1+im2;
    data[neg+1]=re2-im1;

  } /* end for (i=1;i<n/2;i++) */ 

  /* set zero frequency terms to zero */
  data[0]=0.0;
  data[1]=0.0;

  /* set nyquist frequency terms to zero */
  data[n]=0.0;
  data[n+1]=0.0;

  /* note: in general, the zero and nyquist freq terms would be:  */
  /* data[0]=re1; data[1]=re2; data[n]=re1; data[n+1]=re2 */


  /* fft fourier components into the time domain and extract data */
  four1(data-1,n,-1);

  for (i=0;i<n;i++) {
    out1[i]=data[2*i];
    out2[i]=data[2*i+1];
  }

  return;
}
/* Function name: simulate_noise() */

void simulate_noise(int n,float delta_t,double *power,double *whiten_out,
		    float *out,int *pseed)
{

  int    i;
  float  prop,re,im;
  double re_w,im_w,real,imag;

  static int last_n=0;
  static float *data;

  /* (re)allocate memory if current n differs from previous n */

  if (n!=last_n) {
    data=(float *)realloc(data,n*sizeof(float));
    if (data==NULL) {
      GR_start_error("simulate_noise()",rcsid,__FILE__,__LINE__);
      GR_report_error("Failed to allocate %d floats.\n",n);
	GR_end_error();
      abort();
    }

    /* save value of n */
    last_n=n;
  }

  /* simulate noise */

  /* positive frequencies */
  for (i=1;i<n/2;i++) {

    /* random variables */
    prop=2.0/((float)n); /* NOTE: factor of 2 needed for realft() */
    prop*=sqrt((double)n/delta_t);
    prop*=sqrt(power[i]/2.0);
    prop/=sqrt(2.0);
    re=prop*gasdev(pseed);
    im=prop*gasdev(pseed);

    /* whiten data */
    re_w=whiten_out[2*i];
    im_w=whiten_out[2*i+1];

    real=re*re_w-im*im_w;
    imag=re*im_w+im*re_w;
    re=real;
    im=imag;

    /* pack fourier components into data array */
    data[2*i]=re;
    data[2*i+1]=im;

  } /* end for (i=1;i<n/2;i++) */ 

  /* set zero frequency term to zero */
  data[0]=0.0;

  /* set nyquist frequency term to zero */
  data[1]=0.0;

  /* fft fourier components into time domain and extract the data */
  realft(data-1,n,-1);
  for (i=0;i<n;i++) out[i]=data[i];

  return;
}
/* Function Name: combine_data() */

#define MAX_SIZE 16

void combine_data(int which,int n,float *in1,float *in2,float *out)
{
  int    i;
  float  s,c;

  static int last_n[MAX_SIZE];  /* NOTE: static integers are guaranteed */
                                /* to be initialized to zero. */
  static float *buf[MAX_SIZE];


  /* print error message if illegal value for which */

  if (which<1 || which>MAX_SIZE) {
    GR_start_error("combine_data()",rcsid,__FILE__,__LINE__);
    GR_report_error("Illegal value of the argument which(=%d). ",which);
	GR_end_error();
    abort();
  }
  else { 
    /* decrement which so that it can be used as an index to the arrays */
    which=which-1;
  }

  if (n!=last_n[which]) {

    /* (re)allocate memory if current n differs from previous n */

    buf[which]=(float *)realloc(buf[which],(3*n/2)*sizeof(float));
    if (buf[which]==NULL) {
      GR_start_error("combine_data()",rcsid,__FILE__,__LINE__);
      GR_report_error("Failed to allocate %d floats.\n",3*n/2);
	GR_end_error();
      abort();
    }

    /* clear buffer */
    for (i=0;i<3*n/2;i++) buf[which][i]=0.0;

    /* save value of n */
    last_n[which]=n;

  } /* end if (n!=last_n[which]) */

  /* combine the two input arrays (offset by n/2) with one */
  /* another and with data leftover in the buffer. */
  for (i=0;i<n/2;i++) {
    s=sin(i*M_PI/((double)n));
    c=cos(i*M_PI/((double)n));

    buf[which][i]+=in1[i]*s;
    buf[which][i+n/2]+=in1[i+n/2]*c+in2[i]*s;
    buf[which][i+n]+=in2[i+n/2]*c;
  }

  /* output the data */
  for (i=0;i<n;i++) out[i]=buf[which][i];

  /* copy data from last 1/3 of buffer to first 1/3 and clear the rest. */
  /* saved for next call. */

  for (i=0;i<n/2;i++) {
    buf[which][i]=buf[which][i+n];
    buf[which][i+n/2]=0.0;
    buf[which][i+n]=0.0;
  }

  return;
}
