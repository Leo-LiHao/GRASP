/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"
static char *rcsid="$Id: signal_analysis.c,v 1.7 1998/01/23 17:59:49 ballen Exp $\
\n$Name: RELEASE_1_9_8 $";
                        
void four1(float *,int,int);
void realft(float *,int,int);

/* Function name: test_data12() */

#define MAX_DATA   32765  /* magnitude of largest allowed data value */

int test_data12(int n,float *data1,float *data2)
{
  int i,test1=1,test2=1;
  float value;

  static int last_n=0;
  static short *datas;

  /* (re)allocate memory if current n differs from previous value */

  if (n!=last_n) {
    datas=(short *)realloc(datas,n*sizeof(short));
    if (datas==NULL) {
      GR_start_error("test_data12()",rcsid,__FILE__,__LINE__);
      GR_report_error("Failed to allocate %d shorts.\n",n);
	GR_end_error();
      abort();
    }

    /* save value of n */
    last_n=n;
  }
   
  /* test data (site 1) */	

  for (i=0;i<n;i++) {
    value=data1[i];
    if (fabs(value)>MAX_DATA) {
      GR_start_error("test_data12()",rcsid,__FILE__,__LINE__);
      GR_report_error("Warning: Data value %e too large for 16 bits!\n",value);
	GR_end_error();
      value=(value>0.0)?MAX_DATA:-MAX_DATA;
    }
    datas[i]=(short)floor(value+0.5);
  }
  test1=is_gaussian(datas,n,-MAX_DATA-1,MAX_DATA+1,0);

  /* test data (site 2) */

  for (i=0;i<n;i++) {
    value=data2[i];
    if (fabs(value)>MAX_DATA) {
      GR_start_error("test_data12()",rcsid,__FILE__,__LINE__);
      GR_report_error("Warning: Data value %e too large for 16 bits!\n",value);
	GR_end_error();
      value=(value>0.0)?MAX_DATA:-MAX_DATA;
    }
    datas[i]=(short)floor(value+0.5);
  }
  test2=is_gaussian(datas,n,-MAX_DATA-1,MAX_DATA+1,0);

  /* return 1 if both sets passed test. */
  /* otherwise, return 0 and display bad set */

  if (test1 && test2) return 1;
  else {
    if (test1==0) printf("Data segment 1 failed Gaussian test!\n");
    if (test2==0) printf("Data segment 2 failed Gaussian test!\n");
    printf("\n");
    return 0;
  }

}
/* Function Name: analyze() */

void analyze(int average,float *in1,float *in2,int n,float delta_t,
	     float omega_0,float f_low,float f_high,double *gamma12,
	     double *whiten1,double *whiten2,
	     int real_time_noise1,int real_time_noise2,
	     double *power1,double *power2,
	     double *signal,double *variance)
{
  int i;
  float delta_f;

  static int last_n=0;
  static double *filter12;
  static double *signal12;

  /* (re)allocate memory if current n differs from previous value */

  if (n!=last_n) {
    filter12=(double *)realloc(filter12,(n/2)*sizeof(double));
    if (filter12==NULL) {
      GR_start_error("analyze()",rcsid,__FILE__,__LINE__);
      GR_report_error("Error: Failed to allocate %d doubles.\n",n/2);
	GR_end_error();
      abort();
    }  
    signal12=(double *)realloc(signal12,(n/2)*sizeof(double));
    if (signal12==NULL) {
      GR_start_error("analyze()",rcsid,__FILE__,__LINE__);
      GR_report_error("Error: Failed to allocate %d doubles.\n",n/2);
	GR_end_error();
      abort();
    }  

    /* save value of n */
    last_n=n;
  }
   
  /* relation between delta_f and delta_t */
  delta_f=(float)(1.0/(n*delta_t));

  /* extract real-time noise power spectra, if desired */
  if (real_time_noise1==1) 
    extract_noise(average,1,in1,n,delta_t,whiten1,power1);
  if (real_time_noise2==1) 
    extract_noise(average,2,in2,n,delta_t,whiten2,power2);

  /* extract cross-correlation spectrum */
  extract_signal(average,in1,in2,n,delta_t,whiten1,whiten2,signal12);

  /* construct optimal filter function  */
  optimal_filter(n/2,delta_f,f_low,f_high,gamma12,power1,power2,filter12);

  /* calculate cross-correlation signal value corresponding to one */
  /* observation period */
  *signal=0.0;
  for (i=0;i<n/2;i++) {
    if (filter12[i]!=0.0) {
      *signal+=signal12[i]*filter12[i]/(n*delta_t);
    }
  }

  /* calculate theoretical variance for the given noise power spectra */
  *variance=calculate_var(n/2,delta_f,omega_0,f_low,f_high,n*delta_t,
			 gamma12,power1,power2);

  return;
}
/* Function Name: extract_noise() */

#define MAX_SIZE 16
 
void extract_noise(int average,int which,float *in,int n,float delta_t,
		   double *whiten_out,double *power)
{
  int i,j;
  float re,im;
  double re_w,im_w,mag;
  double real,imag;

  static int last_n[MAX_SIZE]; /* NOTE: static integers are guaranteed */
                                /* to be initialized to zero. */
  static float *buf[MAX_SIZE];
  static float *window,*data;
  static double *np[2];


  /* print error message if illegal value for which */

  if (which<1 || which>MAX_SIZE) {
    GR_start_error("extract_noise()",rcsid,__FILE__,__LINE__);
    GR_report_error("Error: Illegal value of the argument which(=%d).\n",which);
	GR_end_error();
    abort();
  }
  else { 
    /* decrement which so that it can be used as an index for the arrays */
    which=which-1;
  }


  if (n!=last_n[which]) {

    /* (re)allocate memory if current n differs from previous value */

    buf[which]=(float *)realloc(buf[which],(3*n/2)*sizeof(float));
    if (buf[which]==NULL) {
      GR_start_error("extract_noise()",rcsid,__FILE__,__LINE__);
      GR_report_error("Failed to allocate %d floats\n",3*n/2);
	GR_end_error();
      abort();
    }
    window=(float *)realloc(window,n*sizeof(float));
    if (window==NULL) {
      GR_start_error("extract_noise()",rcsid,__FILE__,__LINE__);
      GR_report_error("Failed to allocate %d floats\n",n);
	GR_end_error();
      abort();
    }
    data=(float *)realloc(data,n*sizeof(float));
    if (data==NULL) {
      GR_start_error("extract_noise()",rcsid,__FILE__,__LINE__);
      GR_report_error("Failed to allocate %d floats\n",n);
	GR_end_error();
      abort();
    }
    np[0]=(double *)realloc(np[0],(n/2)*sizeof(double));
    if (np[0]==NULL) {
      GR_start_error("extract_noise()",rcsid,__FILE__,__LINE__);
      GR_report_error("Failed to allocate %d doubles\n",n/2);
	GR_end_error();
      abort();
    }
    np[1]=(double *)realloc(np[1],(n/2)*sizeof(double));
    if (np[1]==NULL) {
      GR_start_error("extract_noise()",rcsid,__FILE__,__LINE__);
      GR_report_error("Failed to allocate %d doubles\n",n/2);
	GR_end_error();
      abort();
    }

    /* clear buffer */
    for (i=0;i<3*n/2;i++) buf[which][i]=0.0;

    /* construct hann window function */
    /* (sqrt(8/3) is the "window squared and summed" factor.) */
    /* (see page 553 in Numerical Recipes in C) */
    for (i=0;i<n;i++)
      window[i]=sqrt(8.0/3.0)*0.5*(1.0-cos(2*M_PI*i/((double)n))); 

    /* save value of n */
    last_n[which]=n;

  } /* end if (n!=last_n[which]) */
   

  /* store input data in last 2/3 of buffer */
  for (i=0;i<n;i++) buf[which][i+n/2]=in[i];

  /* analyze the data using the overlapped data segment technique */
  /* (j=0,1 correspond to the two overlapped segments) */

  for (j=0;j<2;j++) {

    /* pack windowed data into an array */
    for (i=0;i<n;i++) data[i]=window[i]*buf[which][i];

    /* fft the windowed data into the freq domain */
    realft(data-1,n,1);

    /* step through the array calculating the noise power spectrum */
    /* at each discrete freq value */

    for (i=0;i<n/2;i++) {

      if (i==0) { /* zero frequency */
	re=data[0];
	im=0.0;
      }
      else { /* positive frequencies */
	re=data[2*i];
	im=data[2*i+1];
      }
      /* NOTE: ignore the (real) nyquist critical frequency component */
      /* contained in data[1]. */

      /* unwhiten the data */
      re_w=whiten_out[2*i];
      im_w=whiten_out[2*i+1];
      mag=re_w*re_w+im_w*im_w;

      real=(re*re_w+im*im_w)/mag;
      imag=(-re*im_w+im*re_w)/mag;
      re=real;
      im=imag;

      /* calculate the value of the noise power spectrum */
      np[j][i]=2.0*(re*re+im*im)*delta_t/n;

    } /* end for (i=0;i<n/2;i++) */

    /* copy last 2/3 of buffers into first 2/3 to overlap the data */
    for (i=0;i<n/2;i++) {
      buf[which][i]=buf[which][i+n/2];
      buf[which][i+n/2]=buf[which][i+n];
    }

  } /* end for (j=0;j<2;j++) */

  /* average the values of the noise power spectrum, if desired */
  for (i=0;i<n/2;i++) {
    if (average==1) power[i]=0.5*(np[0][i]+np[1][i]);
    else power[i]=np[1][i];
  }

  return;
}
/* Function Name: extract_signal() */
              
void extract_signal(int average,float *in1,float *in2,int n,float delta_t,
		    double *whiten1,double *whiten2,double *signal12)
{
  int i,j;
  int pos,neg;
  float gpr,gpi,gnr,gni;
  float re1,im1,re2,im2;
  double re_w1,im_w1,re_w2,im_w2,mag1,mag2;
  double real,imag;

  static int last_n=0;
  static float *buf1,*buf2,*window,*data;
  static double *sig12[2];


  if (n!=last_n) {

    /* (re)allocate memory if current n differs from previous value */

    buf1=(float *)realloc(buf1,(3*n/2)*sizeof(float));
    if (buf1==NULL) {
      GR_start_error("extract_signal()",rcsid,__FILE__,__LINE__);
      GR_report_error("Failed to allocate %d floats.\n",3*n/2);
	GR_end_error();
      abort();
    }
    buf2=(float *)realloc(buf2,(3*n/2)*sizeof(float));
    if (buf2==NULL) {
      GR_start_error("extract_signal()",rcsid,__FILE__,__LINE__);
      GR_report_error("Failed to allocate %d floats.\n",3*n/2);
	GR_end_error();
      abort();
    }
    window=(float *)realloc(window,n*sizeof(float));
    if (window==NULL) {
      GR_start_error("extract_signal()",rcsid,__FILE__,__LINE__);
      GR_report_error("Failed to allocate %d floats.\n",n);
	GR_end_error();
      abort();
    }
    data=(float *)realloc(data,2*n*sizeof(float));
    if (data==NULL) {
      GR_start_error("extract_signal()",rcsid,__FILE__,__LINE__);
      GR_report_error("Failed to allocate %d floats.\n",2*n);
	GR_end_error();
      abort();
    }
    sig12[0]=(double *)realloc(sig12[0],(n/2)*sizeof(double));
    if (sig12[0]==NULL) {
      GR_start_error("extract_signal()",rcsid,__FILE__,__LINE__);
      GR_report_error("Failed to allocate %d doubles.\n",n/2);
	GR_end_error();
     abort();
    }
    sig12[1]=(double *)realloc(sig12[1],(n/2)*sizeof(double));
    if (sig12[1]==NULL) {
      GR_start_error("extract_signal()",rcsid,__FILE__,__LINE__);
      GR_report_error("Failed to allocate %d doubles.\n",n/2);
	GR_end_error();
      abort();
    }

    /* clear buffers */
    for (i=0;i<3*n/2;i++) buf1[i]=buf2[i]=0.0;

    /* construct hann window function */
    /* (sqrt(8/3) is the "window squared and summed" factor.) */
    /* (see page 553 in Numerical Recipes in C) */
    for (i=0;i<n;i++)
      window[i]=sqrt(8.0/3.0)*0.5*(1.0-cos(2*M_PI*i/((double)n))); 

    /* save value of n */
    last_n=n;

  } /* end if (n!=last_n) */
   

  /* store input data in last 2/3 of buffers */
  for (i=0;i<n;i++) {
    buf1[i+n/2]=in1[i];
    buf2[i+n/2]=in2[i];
  }

  /* analyze the data using the overlapped data segment technique */
  /* (j=0,1 correspond to the two overlapped segments) */

  for (j=0;j<2;j++) {

    /* pack windowed data into rea/imag parts of a single array */
    for (i=0;i<n;i++) {
      data[2*i]=window[i]*buf1[i];
      data[2*i+1]=window[i]*buf2[i];
    }

    /* fft the windowed data into the freq domain */
    four1(data-1,n,1);

    /* step through the array calculating the spectrum of the */
    /* cross-correlation signal at each discrete freq value */

    for (i=0;i<n/2;i++) {

      if (i==0) { /* zero frequency */
	re1=data[0];
	re2=data[1];
	im1=im2=0.0;
      }
      else { /* positive frequencies */
	pos=2*i;
	neg=2*(n-i);

	gpr=data[pos];
	gpi=data[pos+1];
	gnr=data[neg];
	gni=data[neg+1];

	re1=0.5*(gpr+gnr);
	im1=0.5*(gpi-gni);
	re2=0.5*(gpi+gni);
	im2=0.5*(gnr-gpr);
      }

      /* unwhiten the data */
      re_w1=whiten1[2*i];
      im_w1=whiten1[2*i+1];
      mag1=re_w1*re_w1+im_w1*im_w1;

      real=(re1*re_w1+im1*im_w1)/mag1;
      imag=(-re1*im_w1+im1*re_w1)/mag1;
      re1=real;
      im1=imag;

      re_w2=whiten2[2*i];
      im_w2=whiten2[2*i+1];
      mag2=re_w2*re_w2+im_w2*im_w2;

      real=(re2*re_w2+im2*im_w2)/mag2;
      imag=(-re2*im_w2+im2*re_w2)/mag2;
      re2=real;
      im2=imag;

      /* calculate the value of the cross-correlation spectrum */
      sig12[j][i]=2.0*(re1*re2+im1*im2)*delta_t*delta_t;

    } /* end for (i=0;i<n/2;i++) */

    /* copy last 2/3 of buffers into first 2/3 to overlap the data */
    for (i=0;i<n/2;i++) {
      buf1[i]=buf1[i+n/2];
      buf1[i+n/2]=buf1[i+n];
      buf2[i]=buf2[i+n/2];
      buf2[i+n/2]=buf2[i+n];
    }

  } /* end for (j=0;j<2;j++) */

  /* average the values of the cross-correlation spectrum, if desired */
  for (i=0;i<n/2;i++) {
    if (average==1) signal12[i]=0.5*(sig12[0][i]+sig12[1][i]);
    else signal12[i]=sig12[1][i];
  }

  return;
}
/* Function Name: optimal_filter() */
                        
void optimal_filter(int n,float delta_f,float f_low,float f_high,
		    double *gamma12,double *power1,double *power2,
		    double *filter12)
{
  int i;
  float  f,f3;
  double sum,norm;


  /* construct the optimal filter function */
  /* (proportional to gamma12(f)/(f^3 P_1(f) P_2(f)) ) */
    
  /* zero frequency */
  filter12[0]=0;

  /* positive frequencies */
  for (i=1;i<n;i++) {
    f=i*delta_f;
    if (f<f_low || f>f_high) filter12[i]=0;
    else {
      f3=f*f*f;
      filter12[i]=gamma12[i]/(f3*power1[i]*power2[i]);
    }
  }

  /* normalize the optimal filter so that <S>=T*Omega_0 */
  sum=0.0;
  for (i=0;i<n;i++) {
    if (filter12[i]!=0.0) sum+=filter12[i]*filter12[i]*power1[i]*power2[i];
  }
  norm=(double)(1.0/delta_f);
  norm/=0.3*HUBBLE*HUBBLE/(M_PI*M_PI);
  norm/=sum;

  for (i=0;i<n;i++) {
    if (filter12[i]!=0.0) filter12[i]*=norm;
  }

  return;
}
/* Function name: calculate_var() */

double calculate_var(int n,float delta_f,float omega_0,
		     float f_low,float f_high,float t,double *gamma12,
		     double *power1,double *power2)
{
  int i;
  float f;
  double factor,f3,f6,f9,f12,p1,p2,g2;
  double int1,int2,int3,int4;
  double variance;

  /* approximate integrals used to calculate the variance */
  int1=int2=int3=int4=0.0;

  for (i=1;i<n;i++) {	/* ignore i=0 contribution to avoid possible */
			/* division by 0 (e.g.,if f_low=0) */
    f=i*delta_f;
    if (f_low<=f && f<=f_high) {
      f3=f*f*f;
      f6=f3*f3;
      f9=f6*f3;
      f12=f6*f6;
      g2=gamma12[i]*gamma12[i];
      p1=power1[i];
      p2=power2[i];

      int1+=delta_f*g2/(f6*p1*p2);
      int2+=delta_f*g2/(f9*p1*p1*p2);
      int3+=delta_f*g2/(f9*p1*p2*p2);
      int4+=delta_f*g2*(1.0+g2)/(f12*p1*p1*p2*p2);

    }
  }

  /* calculate the theoretical variance */
  factor=10.0*M_PI*M_PI/(3.0*HUBBLE*HUBBLE);
  variance=0.5*t*(factor*factor*int1+
                  omega_0*factor*int2+
                  omega_0*factor*int3+
                  omega_0*omega_0*int4)/(int1*int1);
 
  return variance;
}
/* Function name: prelim_stats() */

#define SNR (1.65) /* min SNR required for detection with 95% confidence */

void prelim_stats(float omega_0,float t,double signal,double variance)
{
  double weight;
  double expt_mean,expt_variance,expt_stddev,expt_snr=0.0;
  double theor_mean,theor_variance,theor_stddev,theor_snr;
  double theor_omega_min_95;

  static int runs_completed=0;
  static double weighted_sum1=0.0,weighted_sum2=0.0;
  static double sum_of_weights=0.0;

  /* increment run_completed whenever prelim_stats() is called */
  runs_completed++;

  /* keep track of weighted mean value, variance, and stddev */
  weight=1.0/variance;
  sum_of_weights+=weight;
  weighted_sum1+=weight*signal;  
  weighted_sum2+=weight*signal*signal;

  /* experimental values */
  expt_mean=weighted_sum1/sum_of_weights;
  expt_variance=weighted_sum2/sum_of_weights-expt_mean*expt_mean;
  expt_stddev=sqrt(expt_variance);
  if (runs_completed>1) {
    expt_snr=sqrt((double)runs_completed)*expt_mean/expt_stddev;
  }

  /* theoretical values */
  theor_mean=omega_0*t;
  theor_variance=runs_completed/sum_of_weights;
  theor_stddev=sqrt(theor_variance);
  theor_snr=sqrt((double)runs_completed)*theor_mean/theor_stddev;
  theor_omega_min_95=SNR*theor_stddev/(sqrt((double)runs_completed)*t);

  /* display the latest values */
  printf("total number of runs completed=%d\n",runs_completed);
  printf("total observation time =%e seconds\n",runs_completed*t);
  printf("signal value=%e\n",signal);
  printf("experimental mean=%e\n",expt_mean);
  printf("experimental stddev=%e\n",expt_stddev);
  printf("experimental SNR=%e\n",expt_snr);
  printf("theoretical mean=%e\n",theor_mean);
  printf("theoretical stddev=%e\n",theor_stddev);
  printf("theoretical SNR=%e\n",theor_snr);
  printf("relative error in SNR=%d percent\n",
	 (int)(100*fabs((theor_snr-expt_snr)/theor_snr)));
  printf("experimental omega_0=%e\n",expt_mean/t);
  printf("theoretical omega_0=%e\n",omega_0);
  printf("theoretical omega_0 for detection with 95 percent confidence=%e\n",
	 theor_omega_min_95);
  printf("\n");

  return;
}
/* Function Name: statistics() */

#define GAUSS 8189

void statistics(float *input,int n,int num_bins)
{
  int i,bin,point;
  float min,max,test,width,area,maxprob;
  float x,y,norm;
  double total,total2;
  double mean,variance,stddev,stderror,snr;
  float *histogram,*gaussian;

  FILE *fp1;
  FILE *fp2;

  fp1=fopen("histogram.dat","w");
  fp2=fopen("gaussian.dat","w");

  /* check input data length */
  if (n<2) {
    GR_start_error("statistics()",rcsid,__FILE__,__LINE__);
    GR_report_error("Length of input data=%d < 2.\n",n);
	GR_end_error();
    abort();
  }

  /* allocate memory */

  histogram=(float*)malloc((2*num_bins)*sizeof(float));
  gaussian=(float*)malloc((2*(GAUSS+3))*sizeof(float));

  /* calculate mean, stddev, and confidence intervals */

  /* work out min and max extents */
  min=1.e38;
  max= -min;
  total=total2=0.0;
  for (i=0;i<n;i++) {
    test=input[i];
    if (test<min)
      min=test;
    if (test>max)
      max=test;
    total+=test;
    total2+=test*test;
  }
  mean=total/n;
  variance=(total2/n-mean*mean);
  stddev=sqrt(variance);
  stderror=stddev/sqrt((double)n);
  snr=mean/stderror;

  /* display values */
  fprintf(stdout,"STATISTICS:\n");
  fprintf(stdout,"Number=%d\n",n);
  fprintf(stdout,"Mean=%e\n",mean);
  fprintf(stdout,"Stddev=%e\n",stddev);
  fprintf(stdout,"S/N ratio=%e\n",snr);
  fprintf(stdout,"Con68l=%e\n",mean-stderror);
  fprintf(stdout,"Con68u=%e\n",mean+stderror);
  fprintf(stdout,"Con90l=%e\n",mean-1.65*stderror);
  fprintf(stdout,"Con90u=%e\n",mean+1.65*stderror);
  fprintf(stdout,"Con95l=%e\n",mean-2.0*stderror);
  fprintf(stdout,"Con95u=%e\n",mean+2.0*stderror);
  fprintf(stdout,"Con99l=%e\n",mean-3.0*stderror);
  fprintf(stdout,"Con99u=%e\n",mean+3.0*stderror);
  fprintf(stdout,"\n");


  /* bin the data */

  /* set up the bins */
  width=(max-min)/num_bins;
  for (i=0;i<num_bins;i++) {
    /* label the cbin ranges with the center value */
    histogram[2*i]=min+width*(i+0.5);

    /* and clear the accumulants */
    histogram[2*i+1]=0.0;
  }

  /* do the binning... */
  for (i=0;i<n;i++) {
    test=input[i];
    bin=num_bins*(test-min)/(max-min);
    if (bin>=num_bins)
      bin=num_bins-1;
    histogram[2*bin+1]++;
  }

  /* normalize distribution to unit area */
  maxprob= -1.e38;
  area=0.0;
  for (i=0;i<num_bins;i++)
    area+=width*histogram[2*i+1]; 
  for (i=0;i<num_bins;i++) {
    histogram[2*i+1]/=area;
    if (histogram[2*i+1]>maxprob)
      maxprob=histogram[2*i+1];
  }

  /* GAUSSIAN DISTRIBUTION (and markers) */

  point=0;
  width=8.0*stddev/GAUSS;
  for (i=0;i<GAUSS;i++) {
    x=mean-4.0*stddev+i*width;

    /* mark one standard deviation below mean */
    if (x>(mean-stddev) && point==0) {
      gaussian[2*i+2*point]=x-0.5*width;
      gaussian[2*i+1+2*point]=0.0;
      point++;
    }

    /* mark mean */
    if (x>mean && point==1) {
      gaussian[2*i+2*point]=x-0.5*width;
      gaussian[2*i+1+2*point]=0.0;
      point++;
    }

    /* mark one standard deviation above mean */
    if (x>(mean+stddev) && point==2) {
      gaussian[2*i+2*point]=x-0.5*width;
      gaussian[2*i+1+2*point]=0.0;
      point++;
    }

    /* put in gaussian probablity distribution */
    y=(x-mean)/stddev;
    norm=1.0/(stddev*sqrt(2.0*M_PI));
    gaussian[2*i+2*point]=x;
    gaussian[2*i+1+2*point]=norm*exp(-0.5*y*y);
  }

  /* write histogram and gaussian distribution to files */

  for (i=0;i<num_bins;i++) {
    fprintf(fp1,"%e %e\n",histogram[2*i],histogram[2*i+1]); 
  }

  for (i=0;i<GAUSS+3;i++) {
    fprintf(fp2,"%e %e\n",gaussian[2*i],gaussian[2*i+1]); 
  }

  /* close files */
  fclose(fp1);
  fclose(fp2);

  return;
}

/* NOTE: THIS SHOULD BE REMOVED EVENTUALLY */
/*
void testme2() {
  float *x;
  int y,z;

  graph(x,y,z);
  return;
}
*/
