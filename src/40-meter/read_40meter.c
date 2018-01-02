/* GRASP: Copyright 1997,1998  Bruce Allen */
static char *rcsid="$Id: read_40meter.c,v 1.14 1998/09/04 18:50:05 ballen Exp $\
\n$Name: RELEASE_1_9_8 $";

#include "grasp.h"

#define BLOCK 1024
#define WORD 2
#define HIRATE 9868.4208984375
#define LORATE 986.8421020507812
#define LBUFFER 10000
#define LOCKL 1
#define LOCKH 10

float sample_rate;

/* min and max macros, from NR appendix B */
static int imaxarg1,imaxarg2;
static int iminarg1,iminarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))


/* check that number of points is a power of two... 
	i=1;j=0;
	while (i!=npoint) {
		i*=2;
		j++;
		if (j>31) {
			GR_start_error("normalize_gw()",rcsid,__FILE__,__LINE__);
			GR_report_error("npoint = %d",npoint);
			GR_report_error("should be an integer power of two\n");
			GR_end_error();
			abort();
		}
	}
*/


void normalize_gw(	FILE *fpsine,	/* name of file to get calibration info from.           */
			int npoint,	/* number of points in time domain sample (power of 2)  */
			float srate,	/* sample rate in Hz of the gw signal stream            */
			float *response /* on return, contains complex response function        */
) {

	float freq,factor,f2,mod2,ssreal,ssimag;
	int j,ir,ii;
	
	/* get interpolated swept sine SS curve (store in response[] for now) */
	calibrate(fpsine,npoint,response,srate,2,0);

	/* assemble all the factors -- square root bandwidth */
	factor=sqrt(9.35);

	/* Bob Spero's constant k */
	factor/=21399.0;

	/* volts IFO/ADC bit */
	factor*=(10.0/2048.0);

	/* from two time derivatives */
	factor/=(-4.0*M_PI*M_PI);
	
	/* set DC response to zero */
	response[0]=response[1]=0.0;

	/* set Nyquist response to zero */
	response[npoint]=response[npoint+1]=0.0;

	/* loop over all frequencies except DC (j=0) & Nyquist (j=npoint/2) */
	for (j=1;j<npoint/2;j++) {
		/* subscripts of real, imaginary parts */
		ir=2*j;
		ii=ir+1;

		/* frequency and frequency^2 */
		freq=srate*j/npoint;
		f2=freq*freq;

		/* real, imaginary parts of swept sine SS */
		ssreal=response[ir];
		ssimag=response[ii];

		/* squared modulus of swept sine SS */
		mod2=ssreal*ssreal+ssimag*ssimag;

		/* response function (one over complex conjugate of SS) */
		response[ir]=factor*ssreal/(f2*mod2);
		response[ii]=factor*ssimag/(f2*mod2);
	}
	return;
}



/* 
   method=0 => rational function interpolation
   method=1 => polynomial interpolation
   order in either case is the order of interpolation
*/

void calibrate(FILE *fpsine,int num,float *complex,float srate,int method,int order) {
 	void hunt(float xx[], unsigned long n, float x, unsigned long *jlo);
 	void polint(float xa[], float ya[], int n, float x, float *y, float *dy);
	void ratint(float xa[], float ya[], int n, float x, float *y, float *dy);
	void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);
	void splint(float xa[], float ya[], float y2a[], int n, float x, float *y);

	float *ssf,*ssr,*ssi,freq,interp_real,interp_imag,*ssr2d=NULL,*ssi2d=NULL;
	float error_real,error_imag;
	int i,k,ssn;
	unsigned long jlo=1;

/* read in the swept sine file */
read_sweptsine(fpsine,&ssn,&ssf,&ssr,&ssi);

/* if using spline interpolation, please set up 2nd derivatives */
if (method==2) {
	ssr2d=(float*)malloc(ssn*sizeof(float));
	ssi2d=(float*)malloc(ssn*sizeof(float));
	spline(ssf-1,ssr-1,ssn,2.e30,2.e30,ssr2d-1);
	spline(ssf-1,ssi-1,ssn,2.e30,2.e30,ssi2d-1);
}

/* step through num/2+1 different (non-negative) frequencies */
for (i=0;i<=num/2;i++) {

	/* frequency at which we want response function */
	freq=(i*srate)/num;

	if (method==2) {
		/* use spline interpolation to estimate value */
		splint(ssf-1,ssr-1,ssr2d-1,ssn,freq,&interp_real);	
		splint(ssf-1,ssi-1,ssi2d-1,ssn,freq,&interp_imag);
	}	
	else {
		/* hunt for nearest entries in table (NR routine) -1 for zero offset */
		hunt(ssf-1,ssn,freq,&jlo);

		/* find where to interpolate in table, see NR in C section 3.4 */
		k=IMIN(IMAX(jlo-(order-1)/2,1),ssn+1-order);

		/* now get the interpolated value (-2 is zero offset, twice) */
		if (method==0) {
			/* rational function interpolation */
			ratint(ssf+k-2,ssr+k-2,order,freq,&interp_real,&error_real);
			ratint(ssf+k-2,ssi+k-2,order,freq,&interp_imag,&error_imag);
		}
		else if (method==1) {
			/* polynomial interpolation */
			polint(ssf+k-2,ssr+k-2,order,freq,&interp_real,&error_real);
			polint(ssf+k-2,ssi+k-2,order,freq,&interp_imag,&error_imag);
		}
		else {
			GR_start_error("calibrate()",rcsid,__FILE__,__LINE__);
			GR_report_error("unrecognized interpolation method %d\n",method);
			GR_end_error();
			abort();
		}
	}

	/* store the interpolated values in arrays */
	complex[2*i]=interp_real;
	complex[2*i+1]=interp_imag;
}

/* if using spline interpolation, free the associated memory */
if (method==2) {
	free(ssr2d);
	free(ssi2d);
}

/* free memory from read_sweptsine() */
free(ssf);
free(ssr);
free(ssi);

return;
}



/* blocksize for allocation of calibration arrays */
#define BLOCKSIZE 1024

/* this function opens a swept sine file for reading... */
void read_sweptsine(	FILE *fpsine,	/* name of the swept sine file */
			int *n,		/* number of lines (entries) in file + 1 for DC */
			float **freq,	/* pointer to array of freq values */
			float **real,	/* pointer to array of real parts (in phase) */
			float **imag	/* pointer to array of imaginary parts (90 degrees out of phase) */
) {

int code,i,npoint=0,nread=0;
float *array[3];

/* check file containing swept sine data */
if (fpsine==NULL) {
	GR_start_error("read_calibrate()",rcsid,__FILE__,__LINE__);
	GR_report_error("file pointed to by first argument is null!\n");
	GR_end_error();
	abort();
}

/* create null pointers for freq, real, imag parts */
for (i=0;i<3;i++)
	array[i]=NULL;

/* now loop, getting one additional input line per iteration */
while (1) {
	/* if we need more memory space, then allocate it */
	if (npoint<nread+1) {
		npoint+=BLOCKSIZE;
		for (i=0;i<3;i++) {
			array[i]=(float *)realloc((void *)array[i],sizeof(float)*npoint);
			if (array[i]==NULL) {
				GR_start_error("read_calibrate()",rcsid,__FILE__,__LINE__);
				GR_report_error("failed to allocate %d floats\n",npoint);
				GR_end_error();
				fclose(fpsine);
				abort();
			}
		}
	}

	if (nread==0) {
		/* put in vanishing DC response! */
		array[0][0]=0.00;
		array[1][0]=0.00;
		array[2][0]=0.00;
		code=3;
	}
	else {
		/* read in the data, 3 items per line: freq, real, imag */
		code=fscanf(fpsine,"%e %e %e",array[0]+nread,array[1]+nread,array[2]+nread);
	}

	/* if we are done reading, then exit */
	if (code==EOF) {
		break;
	}

	/* if we did read, but not three items, complain and abort */
	if (code!=3) {
		GR_start_error("read_calibrate()",rcsid,__FILE__,__LINE__);
		GR_report_error("problem on line %d of file pointed to by first argument\n",nread+1);
		GR_end_error();
		fclose(fpsine);
		abort();
	}

	/* increment the number of items read, and continue... */
	nread++;
}

/* If we have read less than 2 lines, something MUST be wrong */
if (nread<2) {
                GR_start_error("read_calibrate()",rcsid,__FILE__,__LINE__);
                GR_report_error("Read less than 2 lines from the swept-sine calibration file!\n");
                GR_report_error("Something must be wrong with the format of this file... exiting.\n");
                GR_end_error();
		exit(1);
}

*n=nread;
*freq=array[0];
*real=array[1];
*imag=array[2];

return;
}


/* This routine gets the next npoints of valid gravity-wave data.  It
   returns 0 if there are no further data streams of the requested
   length.  It returns 1 if starting on a new stream of locked data; it
   returns 2 if the data is part of an on-going continuous locked
   sequence.

   Upon return, *tstart is set to the time stamp of the 0th data
   point.  The data is placed into the location pointed to by location;
   the calling process is assumed to do its own memory management.  The
   integer *rem is the number of points of data remaining in the locked
   sequence, following this call.

   If seek=1, then the routine does not actually copy the data.  This
   is useful for moving a certain distance into the data stream, for
   example.

*/

int get_data(	FILE *fpdata,		/* pointer to data file                   */
		FILE *fplock,		/* name of file containing lock data      */
		float *tstart,		/* returns time of zeroth point           */
		int npoints,		/* number of points desired               */
		short *location,	/* pointer to where data should be put    */
		int *rem,		/* the number of remaining points of data */
		float *fastrate,	/* sample rate of the fast channel, Hz    */
		int seek		/* if 1, seek but don't get data          */
) {
	int mystatus=2,length,eo,settime=1;
	float t1,t2,slowrate;
	static float tret;
	static int so,sb,eb,cb=-1,n,nalloc,remaining=-1,nextn=INT_MAX;
	static int thispage=0;
	static short *data=NULL;
	struct ld_binheader bin_header;
	struct ld_mainheader main_header;

	/* get new locked stretch if needed ... */
	while (remaining<npoints) {
		remaining=find_locked(fplock,&so,&sb,&eo,&eb,&t1,&t2,&slowrate);
		/* we have exhausted the data, please return now... */
		if (remaining==0)
			return 0;
		/* gw data has ten times the sample rate of lock data */
		remaining*=10;
		so*=10;
		mystatus=1;
	}

	/* Number of remaining points after function returns */
	remaining-=npoints;
	*rem=remaining;

	/* move to one block before desired starting point */

	while (cb<(sb-1)) {
		/* now SEEK (no read!) to one block before correct block */
		nextn=read_block(fpdata,&data,&n,&tret,fastrate,1,&nalloc,1,&bin_header,&main_header);
		cb++;
	}
	if (cb<sb)
		thispage=0;

	if (mystatus==1 && cb==sb && thispage)
		thispage=n-so;

 	/* while more points are  needed */
	while (npoints>0) {
		/* if data on this page  */
		if (thispage) {
			/* number of points in data array */
			if (settime) {
				*tstart=tret+so/(*fastrate);
				settime=0;
			}
			length=n-so;
			if (npoints<length)
				length=npoints;
			/* copy data array->output array */
			if (!seek)
				memcpy(location,(char *)(data+so),(size_t)sizeof(short)*length);
			location+=length;
			thispage-=length;
			if (thispage<0) {
				GR_start_error("get_data()",rcsid,__FILE__,__LINE__);
				GR_report_error("trying to read page with ");
				GR_report_error("%d entries (aborting...)\n",thispage);
				GR_end_error();
				abort();
			}
			npoints-=length;
			if (thispage==0)
				so=0;
			else
				so+=length;
		}
		else {
			/* put data onto page */
			if (so!=0 || npoints<nextn) {
				/* copy data into array */
				if (seek && npoints<=nextn)
					nextn=read_block(fpdata,&data,&n,&tret,fastrate,1,&nalloc,0,&bin_header,&main_header);
				else
					nextn=read_block(fpdata,&data,&n,&tret,fastrate,1,&nalloc,seek,&bin_header,&main_header);
				cb++;
				thispage=n-so;
				if (settime) {
					*tstart=tret+so/(*fastrate);
					settime=0;
				}
			}
			else {
				/* or efficiently copy directly to required location... */
				nextn=read_block(fpdata,&location,&n,&tret,fastrate,0,&npoints,seek,&bin_header,&main_header);
				cb++;
				location+=n;
				npoints-=n;
				thispage=0;
				if (settime) {
					*tstart=tret+so/(*fastrate);
					settime=0;
				}
			}
		}

		/* check that current block still nicely in range... */
		if (cb<sb || cb>eb) {
			GR_start_error("get_data()",rcsid,__FILE__,__LINE__);
			GR_report_error("WARNING! ");
			GR_report_error("reading block %d (data starts on block %d ends on %d)\n",cb,sb,eb);
			GR_end_error();
		}
	}
	return mystatus;
}

		
int get_data2(	FILE *fpdata,		/* pointer to data file                   */
		FILE *fplock,		/* name of file containing lock data      */
		float *tstart,		/* returns time of zeroth point           */
		int npoints,		/* number of points desired               */
		short *location,	/* pointer to where data should be put    */
		int *rem,		/* the number of remaining points of data */
		float *fastrate,	/* sample rate of the fast channel, Hz    */
		int seek		/* if 1, seek but don't get data          */
) {
	int mystatus=2,length,eo,settime=1;
	float t1,t2,slowrate;
	static float tret;
	static int so,sb,eb,cb=-1,n,nalloc,remaining=-1,nextn=INT_MAX;
	static int thispage=0;
	static short *data=NULL;
	struct ld_binheader bin_header;
	struct ld_mainheader main_header;

	/* get new locked stretch if needed ... */
	while (remaining<npoints) {
		remaining=find_locked(fplock,&so,&sb,&eo,&eb,&t1,&t2,&slowrate);
		/* we have exhausted the data, please return now... */
		if (remaining==0)
			return 0;
		/* gw data has ten times the sample rate of lock data */
		remaining*=10;
		so*=10;
		mystatus=1;
	}

	/* Number of remaining points after function returns */
	remaining-=npoints;
	*rem=remaining;

	/* move to one block before desired starting point */

	while (cb<(sb-1)) {
		/* now SEEK (no read!) to one block before correct block */
		nextn=read_block(fpdata,&data,&n,&tret,fastrate,1,&nalloc,1,&bin_header,&main_header);
		cb++;
	}
	if (cb<sb)
		thispage=0;

	if (mystatus==1 && cb==sb && thispage)
		thispage=n-so;

 	/* while more points are  needed */
	while (npoints>0) {
		/* if data on this page  */
		if (thispage) {
			/* number of points in data array */
			if (settime) {
				*tstart=tret+so/(*fastrate);
				settime=0;
			}
			length=n-so;
			if (npoints<length)
				length=npoints;
			/* copy data array->output array */
			if (!seek)
				memcpy(location,(char *)(data+so),(size_t)sizeof(short)*length);
			location+=length;
			thispage-=length;
			if (thispage<0) {
				GR_start_error("get_data2()",rcsid,__FILE__,__LINE__);
				GR_report_error("trying to read page with ");
				GR_report_error("%d entries (aborting...)\n",thispage);
				GR_end_error();
				abort();
			}
			npoints-=length;
			if (thispage==0)
				so=0;
			else
				so+=length;
		}
		else {
			/* put data onto page */
			if (so!=0 || npoints<nextn) {
				/* copy data into array */
				if (seek && npoints<=nextn)
					nextn=read_block(fpdata,&data,&n,&tret,fastrate,1,&nalloc,0,&bin_header,&main_header);
				else
					nextn=read_block(fpdata,&data,&n,&tret,fastrate,1,&nalloc,seek,&bin_header,&main_header);
				cb++;
				thispage=n-so;
				if (settime) {
					*tstart=tret+so/(*fastrate);
					settime=0;
				}
			}
			else {
				/* or efficiently copy directly to required location... */
				nextn=read_block(fpdata,&location,&n,&tret,fastrate,0,&npoints,seek,&bin_header,&main_header);
				cb++;
				location+=n;
				npoints-=n;
				thispage=0;
				if (settime) {
					*tstart=tret+so/(*fastrate);
					settime=0;
				}
			}
		}

		/* check that current block still nicely in range... */
		if (cb<sb || cb>eb) {
			GR_start_error("get_data2()",rcsid,__FILE__,__LINE__);
			GR_report_error("WARNING! ");
			GR_report_error("reading block %d (data starts on block %d ends on %d)\n",cb,sb,eb);
			GR_end_error();
		}
	}
	return mystatus;
}

		

		

 /*


 find_locked() returns 0 if it fails to find a locked segment and
 returns n>0 if it does find a locked segment. n is the actual number
 of data points in the locked segment.

 If find_locked() finds a locked segment then *s_block is the absolute
 data block (starting from zero) where the locked segment may be found,
 and *s_offset is the offset number of elements (starting from zero)
 into that data block.  Similarly, *e_block and *e_offset label the end
 of the locked segment *s_offset points the the first locked data
 point, and *e_offset points to the last locked data point.

*/

int find_locked(FILE *fp, int *s_offset,int *s_block,int *e_offset,int *e_block,float *tstart,float *tend,float *srate) {
	int i,code,nsave;
	static float ts,te;
	float tesave;
	static int ls_off,ls_blo,finished;
	static int le_off=-1,le_blo=-1,n=0,nalloc=0;
	static short *data=NULL;
	int totaln;
	struct ld_binheader bin_header;
	struct ld_mainheader main_header;

	/* if no data remains, then */
	if (finished) return 0;

	/* if no lock data, complain */
	if (fp==NULL) {
		GR_start_error("find_locked()",rcsid,__FILE__,__LINE__);
		GR_report_error("data file pointed to by first argument not open!\n");
		GR_end_error();
		abort();
	}
	
	/* We are searching from one after the last end... */
	i=le_off+1;
	ls_blo=le_blo;
	ls_off=-1;
	ts=te;

	while (ls_off<0) {
		/* move pointer along, looking for first locked point */
		while (i<n && (data[i]<LOCKL || data[i]>LOCKH)) i++;
		if (i<n)
			/* we found first locked point */
			ls_off=i;
		else {
			/* read in some new data */
			code=read_block(fp,&data,&n,&ts,srate,1,&nalloc,0,&bin_header,&main_header);
			if (!code) {
				/* no new data, so end gracefully */
				if (nalloc)
					free(data);
				return 0;
			}
			else {
				/* and start looking again for the first locked point */
				ls_blo++;
				i=0;
			}
		}
	}

	/* We've found the first locked point: output location */
	*s_offset=ls_off;
	*s_block=ls_blo;
	*tstart=ts;
	totaln=n;

	/* now look for end, forward from the start... */
	le_blo=ls_blo;
	te=ts;
	i=ls_off;
	le_off=-1;
	while (le_off<0) {
		/* move pointer along, looking for first unlocked point */
		while (i<n && data[i]>=LOCKL && data[i]<=LOCKH)
			i++;
		if (i<n)
			/* we have found an out-of-lock point */
			le_off=i;
		else {
			/* before getting more data, save time, length */
			nsave=n;
			tesave=te;
			code=read_block(fp,&data,&n,&te,srate,1,&nalloc,0,&bin_header,&main_header);
			totaln+=n;
			if (code) {
				le_blo++;
				i=0;
				/* see if we are missing more than 3 data points */
				if (te-tesave>(nsave+3.0)/(double)*srate) {
					GR_start_error("find_locked()",rcsid,__FILE__,__LINE__);
					GR_report_error("Gap in data from t=%f to t=%f\n",
						tesave+nsave/(double)*srate,te);
					GR_end_error();
					/* in this case, output end of last block */
					*e_offset=nsave-1;
					*e_block=le_blo-1;
					*tend=tesave;
					/* next time round, begin search on new data block */
					le_off=-1;
					return (totaln-n-ls_off);
				}
			}
			else {
				/* we are out of data, please exit gracefully... */
				if (nalloc)
					free(data);

				/* set up so that next call returns zero */
				finished=1;

				/* put in offset of final point */
				le_off=n;
			}
		}
	}

	/* output the end-point data */
	le_off--;
	*e_offset=le_off;
	*e_block=le_blo;
	*tend=te;
	return (totaln-ls_off-(n-1-le_off));
}


int find_locked2(FILE *fp, int *s_offset,int *s_block,int *e_offset,int *e_block,float *tstart,float *tend,float *srate) {
	int i,code,nsave;
	static float ts,te;
	float tesave;
	static int ls_off,ls_blo,finished;
	static int le_off=-1,le_blo=-1,n=0,nalloc=0;
	static short *data=NULL;
	int totaln;
	struct ld_binheader bin_header;
	struct ld_mainheader main_header;

	/* if no data remains, then */
	if (finished) return 0;

	/* if no lock data, complain */
	if (fp==NULL) {
		GR_start_error("find_locked2()",rcsid,__FILE__,__LINE__);
		GR_report_error("data file pointed to by first argument not open!\n");
		GR_end_error();
		abort();
	}
	
	/* We are searching from one after the last end... */
	i=le_off+1;
	ls_blo=le_blo;
	ls_off=-1;
	ts=te;

	while (ls_off<0) {
		/* move pointer along, looking for first locked point */
		while (i<n && (data[i]<LOCKL || data[i]>LOCKH)) i++;
		if (i<n)
			/* we found first locked point */
			ls_off=i;
		else {
			/* read in some new data */
			code=read_block(fp,&data,&n,&ts,srate,1,&nalloc,0,&bin_header,&main_header);
			if (!code) {
				/* no new data, so end gracefully */
				if (nalloc)
					free(data);
				return 0;
			}
			else {
				/* and start looking again for the first locked point */
				ls_blo++;
				i=0;
			}
		}
	}

	/* We've found the first locked point: output location */
	*s_offset=ls_off;
	*s_block=ls_blo;
	*tstart=ts;
	totaln=n;

	/* now look for end, forward from the start... */
	le_blo=ls_blo;
	te=ts;
	i=ls_off;
	le_off=-1;
	while (le_off<0) {
		/* move pointer along, looking for first unlocked point */
		while (i<n && data[i]>=LOCKL && data[i]<=LOCKH)
			i++;
		if (i<n)
			/* we have found an out-of-lock point */
			le_off=i;
		else {
			/* before getting more data, save time, length */
			nsave=n;
			tesave=te;
			code=read_block(fp,&data,&n,&te,srate,1,&nalloc,0,&bin_header,&main_header);
			totaln+=n;
			if (code) {
				le_blo++;
				i=0;
				/* see if we are missing more than 3 data points */
				if (te-tesave>(nsave+3.0)/(double)*srate) {
					GR_start_error("find_locked2()",rcsid,__FILE__,__LINE__);
					GR_report_error("Gap in data from t=%f to t=%f\n",
						tesave+nsave/(double)*srate,te);
					/* in this case, output end of last block */
					GR_end_error();
					*e_offset=nsave-1;
					*e_block=le_blo-1;
					*tend=tesave;
					/* next time round, begin search on new data block */
					le_off=-1;
					return (totaln-n-ls_off);
				}
			}
			else {
				/* we are out of data, please exit gracefully... */
				if (nalloc)
					free(data);

				/* set up so that next call returns zero */
				finished=1;

				/* put in offset of final point */
				le_off=n;
			}
		}
	}

	/* output the end-point data */
	le_off--;
	*e_offset=le_off;
	*e_block=le_blo;
	*tend=te;
	return (totaln-ls_off-(n-1-le_off));
}




/* This routine reads in a block of data, including the main header,
   binary header, and data.  It returns 0 on failure (no data
   remaining).  It returns -1 if there will be no further data after
   the current block.  It returns n>0 if the next block contains n data
   points.

   If allocate==0 then the program just puts the data *here, and checks
   to make sure that *nalloc is greater than the number of shorts
   required.  If not, it prints an error and exits.  If allocate!=0
   then the program puts the data *here if *nalloc is larger than the
   storage space required.  Otherwise the program uses realloc to
   create a new data space.  Note that it will properly malloc for the
   first time if *nalloc=0.  However in this case be sure that *here is
   not pointing to a block of empty memory or its address will be lost
   forever and the memory will drift on the heap...
   
*/
int read_block(
	FILE *fp,         /* points to file                             */
	short **here,     /* place to put valid data                    */
	int *n,           /* number of data points                      */
	float *tstart,    /* time at start of block                     */
	float *srate ,    /* sample rate for data in Hz                 */
	int allocate,     /* if 0, don't malloc, if 1 malloc or realloc */
	int *nalloc,       /* number of shorts of memory allocated       */
	int seek,          /* if 1, don't read data, just seek...        */
	struct ld_binheader *bin_header,  /* pointer to bin header     */
	struct ld_mainheader *main_header/* pointer to main header     */

) {
	struct ld_mainheader main_header2;
	int code,localn,ntoswap;
	static int bigendian=-1;
	void swap2(int,char *);
	void swap4(int,char *);
	union {
		short shortint;
		char  twobyte[2];
	} test;

	/* check that the structure sizes are OK */
	if (sizeof(struct ld_mainheader)!=48 || sizeof(struct ld_binheader)!=8) {
		GR_start_error("read_block()",rcsid,__FILE__,__LINE__);
		GR_report_error("The size of the header structures is wrong.\n");
		GR_report_error("Apparently your machine does not use 4 byte ints and 4 byte floats.\n");
		GR_report_error("Change grasp.h to get the header sizes for the old 40-meter data format.\n");
		GR_report_error("You need to have sizeof(struct ld_mainheader)=48 and sizeof(struct ld_binheader)=8.\n");
		GR_report_error("aborting....\n");
		GR_end_error();
		abort();
	}

	/* first pass through, determine the byte ordering of the machine */
	if (bigendian==-1) {
		test.shortint=0x100;
		if (test.twobyte[1]) {
			bigendian=0;
			GR_start_error("read_block()",rcsid,__FILE__,__LINE__);
			GR_report_error("Your hardware platform is NOT big-endian.\n");
			GR_report_error("All data will be byte-swapped before being processed.\n");
			GR_end_error();
		}
		else
			bigendian=1;
	}


	/* Read the main header, return if incomplete */
	code=fread((void*)main_header,(size_t)sizeof(struct ld_mainheader),(size_t)1,fp);
	if (code!=1) return 0;

	/* Read the binary header, return if incomplete */
	code=fread((void*)bin_header,(size_t)sizeof(struct ld_binheader),(size_t)1,fp);
	if (code!=1) return 0;


	if (!bigendian) {
		ntoswap=sizeof(struct ld_mainheader)/4;
		swap4(ntoswap,(char *)main_header);
		ntoswap=sizeof(struct ld_binheader)/4;
		swap4(ntoswap,(char *)bin_header);
	}

	/* reality check... */
	if (bin_header->datarate!=LORATE && bin_header->datarate!=HIRATE) {
		GR_start_error("read_block()",rcsid,__FILE__,__LINE__);
		GR_report_error("Warning - unknown datarate %f\n",
			bin_header->datarate);
		GR_report_error("malformed data? %f\n",bin_header->datarate);
		GR_end_error();
	}

	/* Set number of shorts  */
	localn=main_header->chunksize/sizeof(short);

	/* Allocate/reallocate storage space, if desired */
	if  (allocate) {
		if (*nalloc==0)
			/* make pointer to nowhere before allocating space first time */
			*here=NULL;
		if (*nalloc<localn) {
			/* allocate storage space, please */
			*here=(short*)realloc(*here,localn*sizeof(short));
			*nalloc=localn;
		}
		if (*here==NULL) {
			GR_start_error("read_block()",rcsid,__FILE__,__LINE__);
			GR_report_error("failed to allocate %d shorts\n",localn);
			GR_end_error();
			abort();
		}
	}
	else if (*nalloc<localn && !seek) {
		/* in this case, we are storing data in existing memory, but there's not enough space! */
		GR_start_error("read_block()",rcsid,__FILE__,__LINE__);
		GR_report_error("can not put %d shorts into %d space\n",localn,*nalloc);
		GR_end_error();
		abort();
	}

	if (seek)
		fseek(fp,sizeof(short)*localn+BLOCK-sizeof(struct ld_binheader)-sizeof(struct ld_mainheader),SEEK_CUR);
	else {

		/* Skip remaining zeros */
		fseek(fp,BLOCK-sizeof(struct ld_binheader)-sizeof(struct ld_mainheader),SEEK_CUR);

		/* and read in the data, returning if incomplete */
		code=fread((void*)*here,(size_t)sizeof(short),(size_t)localn,fp);
		if (code!=localn)
			return 0;

		/* swap bytes if needed */
		if (!bigendian)
			swap2(localn,(char *)*here);
	}

	/* Sucess!  Set outputs, tell world about it... */
	*n=localn;
	*tstart=bin_header->elapsed_time;
	*srate=bin_header->datarate;

	/* Get chunksize from next block, return -1 if no further blocks */
	code=fread((void*)&main_header2.chunksize,(size_t)4,(size_t)1,fp);
	if (code!=1)
		return -1;

	/* seek back to start of the block */
	fseek(fp,-4,SEEK_CUR);

	if (!bigendian)
		swap4(1,(char *)&main_header2.chunksize);

	/* get the size of the *next* data block */
	localn=main_header2.chunksize/sizeof(short);
	
	return localn;
}

void swap2(int n,char *here) {

	char *next;
	char temp;
	int i;

	for (i=0;i<n;i++) {
		next=here+1;
		temp=*here;
		*here=*next;
		*next=temp;
		here+=2;
	}
	return;
}

void swap4(int n, char *b1) {
	char *b2,*b3,*b4;
	char temp;
	int i;

	for (i=0;i<n;i++) {
		b2=b1+1;
		b3=b2+1;
		b4=b3+1;

		temp=*b1;
		*b1=*b4;
		*b4=temp;
		temp=*b2;
		*b2=*b3;
		*b3=temp;
		b1+=4;
	}
	return;
}
