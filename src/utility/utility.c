/* GRASP: Copyright 1997,1998  Bruce Allen */
static char *rcsid="$Id: utility.c,v 1.23 1999/07/01 20:51:41 ballen Exp $\
\n$Name: RELEASE_1_9_8 $";

#include "grasp.h"
#ifndef _WIN32 
#include <sys/utsname.h>
#endif 

/* when called with pointers to arrays c,a,b, this routine
   takes sets c = a x b where a,b,c are complex numbers */

void product(float *c,float *a, float *b,int ncomplex) {
	float ar,ai,br,bi;

	while (ncomplex-->0) {

		ar=*a++;
		ai=*a++;
		br=*b++;
		bi=*b++;

		*c++=ar*br-ai*bi;
		*c++=ar*bi+ai*br;
	}
	return;
}

/* when called with pointers to arrays c,a,b, this routine
   sets c = a / b where a,b,c are complex numbers */

void ratio(float *c,float *a, float *b,int ncomplex) {
	float ar,ai,br,bi,mod;

	while (ncomplex-->0) {

		ar=*a++;
		ai=*a++;
		br=*b++;
		bi=*b++;

		mod=1.0/(br*br+bi*bi);

		*c++=mod*(ar*br+ai*bi);
		*c++=mod*(ai*br-ar*bi);
	}
	return;
}



/* when called with pointers to arrays c,a,b, this routine
   sets c = a x b* where a,b,c are complex numbers,
   and * is conjugate. */

void productc(float *c,float *a, float *b,int ncomplex) {
	float ar,ai,br,bi;

	while (ncomplex-->0) {

		ar=*a++;
		ai=*a++;
		br=*b++;
		bi=*b++;

		*c++=ar*br+ai*bi;
		*c++=ai*br-ar*bi;
	}
	return;
}

void avg_spec(float *data, float *average, int npoint, int *reset,
		  float srate, float decaytime, int windowtype, int overlap)
{
    void realft(float*,unsigned long,int);
    static float *saved=NULL,*working=NULL,*window=NULL;
    static double wss,decay,lastnorm;
    float spec,*dataptr;
    double win,delta,norm;
    int i,ipn,im,re;

    /* Keep the data */
    dataptr=data;

    /* if we are resetting, then.. */
    if (*reset==0) {
	/* allocate space for saved data */
	saved=(float *)realloc(saved,sizeof(float)*npoint/2);
	if (saved==NULL) {
	    GR_start_error("avg_spec()",rcsid,__FILE__,__LINE__);
	    GR_report_error("failed to allocate %d floats\n",npoint);
	    GR_end_error();
	    abort();
	}
	
	/* allocate working space */
	working=(float *)realloc(working,sizeof(float)*npoint);
	if (working==NULL) {
	    GR_start_error("avg_spec()",rcsid,__FILE__,__LINE__);
	    GR_report_error("failed to allocate %d floats\n",npoint);
	    GR_end_error();
	    abort();
	}
	
	/* allocate space for window function */
	window=(float *)realloc(window,sizeof(float)*npoint);
	if (window==NULL) {
	    GR_start_error("avg_spec()",rcsid,__FILE__,__LINE__);
	    GR_report_error("failed to allocate %d floats\n",npoint);
	    GR_end_error();
	    abort();
	}
	
	/* set average to zero and put data in saved space */
	for (i=0;i<npoint/2;i++) {
	    saved[i]= *(dataptr++);
	    average[i]=0.0;
	}

	/* compute windows  Note n = 2 x npoint */
	wss=0.0;
	for (i=0;i<npoint;i++) {
	    switch (windowtype) {
		/* rectangular (no) window */
	    case 0: win=1.0;
		break;
		/* Hann window */
	    case 1:	win=0.5*(1.0-cos(2.0*M_PI*i/npoint));
		break;
		/* Welch window */
	    case 2: win=(2.0*i)/npoint-1.0;
		win=1.0-win*win;
		break;
		/* Bartlett window */
	    case 3: win=1.0-fabs((2.0*i)/npoint-1.0);
		break;
	    default: 
		GR_start_error("avg_spec()",rcsid,__FILE__,__LINE__);
		GR_report_error("don't recognize windowtype=%d\n",
				windowtype);
		GR_end_error();
		abort();
		break;
	    }
	    /*  Wss differs by factor of N from Numerical Recipes */
	    wss+=win*win;
	    window[i]=win;
	}
	/* normalization: gives units of data^2/Hz to output, one-sided power spectrum */
	wss=2.0/(wss*srate);
	
	/* set decay constant for averaging, from decay time */
	decay=exp(-1.0*npoint/(srate*decaytime));
	
	/* resets averaging over power spectrum */
	lastnorm=0.0;
	
	/* say that we are now resetting */
	*reset=1;

	/* indicate that the first call has no data to overlap */
	if(overlap)
	    overlap--;

    }else{

	/* If not overlapping,  must fill saved with first half of data */
	if(!(overlap))
	    for (i=0;i<npoint/2;i++) saved[i]= *(dataptr++);

    }

    /* Loop over the data twice if overlapping */
    while (overlap-->=0){
	/* window the data data */
	for (i=0;i<npoint/2;i++) {
	    ipn=i+npoint/2;
	    
	    /* multiply window function times saved data */
	    working[i]=window[i]*saved[i];
	
	    /* multiply window function times new data */
	    working[ipn]=window[ipn]*(*dataptr);
	
	    /* save new data */
	    saved[i]= *(dataptr++);
	}
    
	/* compute FFT of windowed data */
	realft(working-1,npoint,1);
    
	/* for averaging... */
	norm=1.0+decay*lastnorm;
	delta=decay*lastnorm/norm;
    
	/* do DC term (don't want Nyquist) */
	average[0]*=delta;
	average[0]+=0.5*wss*working[0]*working[0]/norm;
    
	/* now step through FFT */
	for (i=1;i<npoint/2;i++) {
	    /* compute subscripts of real, imag parts */
	    im=(re=i+i)+1;
	    
	    /* compute the spectrum */
	    spec=wss*(working[re]*working[re]+working[im]*working[im]);
	    /* do moving average of spectrum */
	
	    average[i]*=delta;
	    average[i]+=spec/norm;
	}
    
	/* set averaging constant for next time */
	lastnorm=norm;
    }
  
    return;
}

/*
This routine looks at the data in the input[0..ninput-1] array.  It
then increments the bins[] array, in histogram style, modifying
bins[0..nbins-1].  output[0] is increased by the number of elements of
input[] in the range [firstbin,firstbin+width), output[1] is increased
by the number of elements of input[] in the range
[firstbin+width,firstbin+2*width) and so on.  *under and *over are
incremented by the number of elements which are under range or over range.
*/
void binner(float *input,int ninput,double *bins,int nbins,
            double firstbin,double width,int *under,int *over){

	int i,offset;

	/* step through the input array */
	for (i=0;i<ninput;i++) {

		/* calculate in which bin the point would belong */
		offset=floor((input[i]-firstbin)/width);

		/* indicate if the bin is below range... */
		if (offset<0)
			(*under)++;
		/* or if the bin is above range ... */
		else if (offset>=nbins)
			(*over)++;
		else
		/* and if it's in range, then bin it! */
			bins[offset]+=1.0;
	}
	return;
}

/*
  Here input points to an array of shorts, of length ninput.  This
  routines takes these values, input[0..ninput-1] and bins them into
  bins.  Each element input[r] increments the value of
  bins[input[r]+offset] by one.
*/
  
	
void binshort(short *input,int ninput,double *bins,int offset){
	int i;
	/* step through the input array */
	for (i=0;i<ninput;i++)  
		/* increment the bin to which the point would belong */
		bins[input[i]+offset]+=1.0;
	return;
}
		
void clear(float *array,int n,int spacing) {
	int i;
	for (i=0;i<n;i++)
		array[i*spacing]=0.0;
	return;
}

void graph(float *array,int n,int spacing) {
	FILE *fp;
	int i;

	/* open a file for writing */
	fp=fopen("temp.graph","w");
        if (fp==NULL) {
	    GR_start_error("graph()",rcsid,__FILE__,__LINE__);
	    GR_report_error("Can't write file \"temp.graph\" in current directory -- sorry!\n");
	    GR_end_error();
	    return;
	}

	/* print data into the file */
	for (i=0;i<n;i++)
		fprintf(fp,"%d\t%e\n",i,array[i*spacing]);

	/* close the file */
	fclose(fp);

	/* start up graphing program with data in the file */
#ifdef __MACOS__
	GnuPlotCommand("plot temp.graph");
#else
	system("xmgr temp.graph 1>/dev/null 2>&1 &");
#endif
	/* and return */
	return;
}

void graph_double(double *array,int n,int spacing) {
	FILE *fp;
	int i;

	/* open a file for writing */
	fp=fopen("temp.graph","w");
        if (fp==NULL) {
	    GR_start_error("graph_double()",rcsid,__FILE__,__LINE__);
	    GR_report_error("Can't write file \"temp.graph\" in current directory -- sorry!\n");
	    GR_end_error();
	    return;
	}

	/* print data into the file */
	for (i=0;i<n;i++)
		fprintf(fp,"%d\t%e\n",i,array[i*spacing]);

	/* close the file */
	fclose(fp);

	/* start up graphing program with data in the file */
#ifdef __MACOS__
	GnuPlotCommand("plot temp.grap");
#else
	system("xmgr temp.graph 1>/dev/null 2>&1 &");
#endif

	/* and return */
	return;
}

void graph_short(short *array,int n) {
	FILE *fp;
	int i;

	/* open a file for writing */
	fp=fopen("temp.graph","w");
        if (fp==NULL) {
	    GR_start_error("graph_short()",rcsid,__FILE__,__LINE__);
	    GR_report_error("Can't write file \"temp.graph\" in current directory -- sorry!\n");
	    GR_end_error();
	    return;
	}

	/* print data into the file */
	for (i=0;i<n;i++)
		fprintf(fp,"%d\t%d\n",i,array[i]);

	/* close the file */
	fclose(fp);

	/* start up graphing program with data in the file */
#ifdef __MACOS__
	GnuPlotCommand("plot temp.grap");
#else
	system("xmgr temp.graph 1>/dev/null 2>&1 &");
#endif

	/* and return */
	return;
}

void sgraph(short *array,int n,char* ext,int filenumber) {
	FILE *fp;
	int i;
	char name[256];

	/* open a file for writing */
	sprintf(name,"%s.%03d",ext,filenumber);
	fp=fopen(name,"w");
        if (fp==NULL) {
	    GR_start_error("sgraph()",rcsid,__FILE__,__LINE__);
	    GR_report_error("Can't write file %s in current directory -- sorry!\n",name);
	    GR_end_error();
	    return;
	}

	/* print data into the file */
	for (i=0;i<n;i++)
		fprintf(fp,"%d\t%d\n",i,array[i]);

	/* close the file */
	fclose(fp);

	/* and return */
	return;
}

void audio(short *array,int n) {
	int i;
	FILE *fp;
	short *temp;
	short store,max=0;
	float norm;

	/* write the data */
	if (NULL==(temp=(short *)malloc(n*sizeof(short)))) {
 		GR_start_error("audio()",rcsid,__FILE__,__LINE__);
		GR_report_error("Unable to allocate %d bytes memory for audio array",n*sizeof(short)); 
 		GR_end_error();
 		return; 
 	}

	/* open a file for writing */
	fp=fopen("temp.au","w");
        if (fp==NULL) {
	    GR_start_error("audio()",rcsid,__FILE__,__LINE__);
	    GR_report_error("Can't write file \"temp.au\" in current directory -- sorry!\n");
	    GR_end_error();
	    return;
	}

	/* magic number */
	i=0x2e736e64;
	fwrite(&i,4,1,fp);

	/* header length */
	i=24;
	fwrite(&i,4,1,fp);

	/* data size (bytes) */
	i=2*n;
	fwrite(&i,4,1,fp);

	/* encoding (16 bit pcm) */
	i=3;
	fwrite(&i,4,1,fp);

	/* sample rate */
	i=9600;
	fwrite(&i,4,1,fp);

	/* number of channels */
	i=1;
	fwrite(&i,4,1,fp);

	/* find max... */
	for (i=0;i<n;i++)
		if (max<(store=abs(array[i])))
			max=store;
	norm = (((float)(SHRT_MAX-1))/(float)max);

	/* write the data */
	for (i=0;i<n;i++) temp[i]=(int)(norm*array[i]);
	fwrite(temp,2,n,fp);
	free(temp);

	fclose(fp);

#ifndef _WIN32
{	/* determine what type of system we are on */
	struct utsname systemtype;
	if (uname(&systemtype)<0) {
 		GR_start_error("audio()",rcsid,__FILE__,__LINE__);
		GR_report_error("Unable to determine host system type: uname() system call failed\n"); 
 		GR_end_error();
 		return; 
 	}

	/* check to see if the operating system is one that we know how to play sound on */
	if (strcmp(systemtype.sysname,"SunOS")==0)
		system("audioplay temp.au");
	else if (strcmp(systemtype.sysname,"Linux")==0)
		system("sox temp.au -U -b -r 8000 -t au /dev/audio");
	else if (strcmp(systemtype.sysname,"IRIX")==0)
		system("sfplay temp.au");
	else {
		/* since we don't recognize the operating system type, print a helpful error message */
		GR_start_error("audio()",rcsid,__FILE__,__LINE__);
		GR_report_error("The audio() function does not know how to play a *.au file\n"); 
		GR_report_error("on the %s operating system, because we don't have a machine\n",systemtype.sysname);
		GR_report_error("of this type available. If you can determine how to do this, please\n"); 
		GR_report_error("email ballen@dirac.phys.uwm.edu with that information, and future\n"); 
		GR_report_error("releases of GRASP will play sound correctly on your OS type.\n"); 
 		GR_end_error();
 	}
}
#else
		GR_start_error("audio()",rcsid,__FILE__,__LINE__);
		GR_report_error("The audio() function does not know how to play a *.au file\n"); 
		GR_report_error("on the WIN32 operating system, because we don't have a machine\n");
		GR_report_error("of this type available. If you can determine how to do this, please\n"); 
		GR_report_error("email ballen@dirac.phys.uwm.edu with that information, and future\n"); 
		GR_report_error("releases of GRASP will play sound correctly on your OS type.\n"); 
 		GR_end_error();
#endif
	return;
}

void sound(short *array,int n,char *ext,int filenumber) {
	int i;
	FILE *fp;
	short *temp;
	short store,max=0;
	float norm;
	char name[256];

	/* write the data */
	if (NULL==(temp=(short *)malloc(n*sizeof(short)))) {
 		GR_start_error("sound()",rcsid,__FILE__,__LINE__);
		GR_report_error("Unable to allocate %d bytes memory for sound array",n*sizeof(short)); 
 		GR_end_error();
 		return; 
 	}

	/* open a file for writing */
	sprintf(name,"%s.%03d.au",ext,filenumber);
	fp=fopen(name,"w");

        if (fp==NULL) {
	    GR_start_error("sound()",rcsid,__FILE__,__LINE__);
	    GR_report_error("Can't write file %s in current directory -- sorry!\n",name);
	    GR_end_error();
	    return;
	}

 	/* magic number */
	i=0x2e736e64;
	fwrite(&i,4,1,fp);

	/* header length */
	i=24;
	fwrite(&i,4,1,fp);

	/* data size (bytes) */
	i=2*n;
	fwrite(&i,4,1,fp);

	/* encoding (16 bit pcm) */
	i=3;
	fwrite(&i,4,1,fp);

	/* sample rate */
	i=9600;
	fwrite(&i,4,1,fp);

	/* number of channels */
	i=1;
	fwrite(&i,4,1,fp);

	/* find max... */
	for (i=0;i<n;i++)
		if (max<(store=abs(array[i])))
			max=store;

	norm = (((float)(SHRT_MAX-1))/(float)max);


	/* write the data */
	for (i=0;i<n;i++) temp[i]=(int)(norm*array[i]);
	fwrite(temp,2,n,fp);
	free(temp);

	fclose(fp);

	return;
}

int is_gaussian(short *array,int n,int min,int max,int print) {

	float erff(float);
	float erffc(float);
	double *histogram,mean,*hcenter;
	float exceed,sum01,sum13,sum35;
	int nbins,retval,i,sigma;
	short value,central_bin;
	static double *storage=NULL;
	static int lastbins=0;

	/* number of bins in histogram */
	nbins=max-min+1;

	/* allocate memory for histogram (factor of 5 for safety!) */
	if (lastbins<nbins && storage!=NULL)
		free(storage);

	if (storage==NULL) {
		storage=(double *)malloc(5*sizeof(double)*nbins);
		lastbins=nbins;
		if (storage==NULL) {
			GR_start_error("is_gaussian()",rcsid,__FILE__,__LINE__);
			GR_report_error("unable to allocate memory ");
			GR_report_error("for histogram of %d distinct values\n",max-min+1);
			GR_end_error();
			abort();
		}
	}


	/* pointer to "zero" bin (use storage+2*lastbins->storage+3*lastbins) */
	histogram=storage+2*lastbins;

	/* clear out the histogram */
	for (i=0;i<nbins;i++) histogram[i]=0.0;


	/* bin the points, and find the mean */
	mean=0.0;
	for (i=0;i<n;i++) {
		value=array[i];
		histogram[value-min]+=1.0;
		mean+=value;
	}
	mean/=n;

	/* and the central bin in which the mean value would fall */
	central_bin=floor(mean+0.5);
	hcenter=histogram+central_bin-min;

	/* robust estimator of sigma, use fact that 68% of points in +- 1 sigma */
	exceed=n*erff(sqrt(0.5));
	sum01=hcenter[0];
	for (i=1;sum01<exceed;i++)
		sum01+=(hcenter[-i]+hcenter[i]);
	sigma=i-1;

	/* now find number of points from one to three sigma  */
	sum13=0.0;
	for (;i<3*sigma;i++)
		sum13+=hcenter[-i]+hcenter[i];

	/* now find number of points from three to five sigma */
	sum35=0.0;
	for (;i<=5*sigma;i++)
		sum35+=hcenter[-i]+hcenter[i];

	/* decide if the distribution is Gaussian.  Number of points in +- 3 sigma  */
	exceed=n*erffc(3.0*sqrt(0.5));
	retval=1;
	if (sum35>exceed+2.0*sqrt(exceed) || sum01+sum13+sum35!=n)
		retval=0;

	/* print info if desired  */
	if (print)
		printf("   Distribution: s= %d, N>3s= %d (expect %d), N>5s= %d (expect 0)\n",
			sigma,(int)sum35,(int)exceed,(int)(n-sum01-sum35-sum13));

	return retval;
}


/* This routine opens files, appending the path determined by the
   environment variable with the character string shortpath.  It tries
   to print intelligent error messages if the environment variable is
   not defined, or if the file pointed to does not exist.
*/

FILE* grasp_open(const char *environment_variable,const char *shortpath,const char *mode) {

	char *pathhead;
	char longpath[256];
	FILE *fp;

	/* get the environment variable giving the data path */
	pathhead=getenv(environment_variable);

	/* if this is null, then print out a warning to the user */
	if (!pathhead) {
		GR_start_error("grasp_open()",rcsid,__FILE__,__LINE__);
		GR_report_error("the environment variable %s must be set.\n",environment_variable);
		GR_report_error("It should point to a directory containing 40-meter data, or parameters.\n");
		GR_report_error("It can be set with a command like:\nsetenv %s /this/is/the/path\n",environment_variable);
		GR_end_error();
		exit(1);
	}

	/* now copy this (head of the data path) into the return */
	strcpy(longpath,pathhead);

	/* and append the remainder of the path */
#ifdef __MACOS__
	strcat(longpath,":");
#elif defined _WIN32
	strcat(longpath,"\\");
#else
	strcat(longpath,"/");
#endif
	strcat(longpath,shortpath);

	/* now open the file */
	fp=fopen(longpath,mode);

	if (fp==NULL) {
		GR_start_error("grasp_open()",rcsid,__FILE__,__LINE__);
		GR_report_error("Unable to open file %s in mode \"%s\"\n",longpath,mode);
		GR_report_error("The current value of environment variable %s is %s\n",environment_variable,pathhead);
		GR_end_error();
		exit(1);
	}

	return fp;
}

	
	
