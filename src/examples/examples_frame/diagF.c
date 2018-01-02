/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"
#define NPOINT 2048        /* determines frequency resolution */
#define SHIFT 1024         /* determines time resolution */
#define CHANNEL 0          /* channel number */        
#define AVG_TIME 10.0      /* seconds */
#define THRESHOLD 6.0      /* for outputting a point */
#define INCLUDE 10.0       /* threshold for including in averaging */
#define SEC_INTO_LOCK 60.0 /* number of seconds into locked segment to start analysis */

int main(int argc,char **argv) {
	void realft(float*,unsigned long,int);
	void clear(float*,int,int);
	FILE *fpout=NULL;
	short *datas;
	int i,code,npoint,im,re,initialize=10,filenumber=0;
	float spec,srate,tstart,*meanp,*meanv,*data,*windata,decaytime;
	float win,freq,deviation,var,*window;
	double decay=0.0,delta,norm,lastnorm=0.0;
	char name[64],*title;
	time_t labtime;
	struct fgetinput fgetinput;
	struct fgetoutput fgetoutput;

	/* number of points to sample and fft (power of 2) */
	npoint=NPOINT;

	/* create arrays */
	datas=(short*)malloc(sizeof(short)*SHIFT);
	data=(float *)malloc(sizeof(float)*npoint);
	windata=(float *)malloc(sizeof(float)*npoint);
	window=(float *)malloc(sizeof(float)*npoint);
	meanp=(float *)malloc(sizeof(float)*npoint/2);
	meanv=(float *)malloc(sizeof(float)*npoint/2);

	/* set up fgetinput for getting needed quantities */
	fgetinput.nchan=1;
	fgetinput.files=framefiles;
	fgetinput.chnames=(char **)malloc(fgetinput.nchan*sizeof(char *));
	fgetinput.locations=(short **)malloc(fgetinput.nchan*sizeof(short *));
	fgetoutput.npoint=(int *)malloc(fgetinput.nchan*sizeof(int));
	fgetoutput.ratios=(int *)malloc(fgetinput.nchan*sizeof(int));
	fgetinput.npoint=SHIFT;
	fgetinput.seek=0;
	fgetinput.calibrate=0;
	fgetinput.locations[0]=datas;

	/* channel names */
	fgetinput.chnames[0]="IFO_DMRO";

	/* set up different cases */
	if (NULL!=getenv("GRASP_REALTIME")) {
		/* 40-meter lab */
		fgetinput.inlock=0;
		fgetinput.chnames[0]=getenv("GRASP_REALTIME");
	}
	else {
		/* Nov 1994 data set */
		fgetinput.inlock=1;
	}

	/* compute Welch window, once */
	for (i=0;i<npoint;i++) {
		win=(2.0*i)/npoint-1.0;
		window[i]=1.0-win*win;
	}

/* This is the main program loop, which aquires data, then filters it */
while ((code=fget_ch(&fgetoutput,&fgetinput))) {
	tstart=fgetoutput.dt;
	srate=fgetoutput.srate;

	/* if just entering a new locked stretch */
	if (code==1) {
		/* decaying exponential average of spectrum: */
		decaytime=AVG_TIME; /* seconds! */
		decay=exp(-1.0*SHIFT/(srate*decaytime));
		initialize=10;

		/* resets averaging over power spectrum */
		lastnorm=0.0;
		clear(meanp,npoint/2,1);
		clear(meanv,npoint/2,1);

		/* clear stored data */
		clear(data,npoint,1);

		if (fpout!=NULL) fclose(fpout);
		sprintf(name,"ch%ddiag.%03d",CHANNEL,filenumber);
		fpout=fopen(name,"w");
 		fprintf(fpout,"@ sets linestyle 0\n");
 		fprintf(fpout,"@ sets symbol 1\n");
 		fprintf(fpout,"@ s%d symbol color %d\n",CHANNEL,CHANNEL+1);
		fprintf(fpout,"@ subtitle \"Time/Frequency statistics for channel %d\"\n",CHANNEL);
		fprintf(fpout,"@ xaxis label \"Time (sec)\"\n");
		fprintf(fpout,"@ yaxis label \"Frequency (Hz)\"\n");
		if (NULL!=getenv("GRASP_REALTIME")) {
			time(&labtime);
			fprintf(fpout,"@ title \"40-meter Lab %.24s\"\n",ctime(&labtime));
		}
		else if ((title=getenv("GRASP_FRAMEPATH")))
			fprintf(fpout,"@ title \"%s\"\n",title);
		filenumber++;
 	}

	/* shift older-time data up to start of array */
	for (i=0;i<npoint-SHIFT;i++)
		data[i]=data[i+SHIFT];

	/* put in the most recent time data */
	for (i=0;i<SHIFT;i++)
		data[i+npoint-SHIFT]=datas[i];

	/* wait a certain amount of time into lock */
	if (fgetoutput.tstart-fgetoutput.lastlock<SEC_INTO_LOCK) continue;

	/* then window it */
	for (i=0;i<npoint;i++)
		windata[i]=window[i]*data[i];

	/* find the FFT */
	realft(windata-1,npoint,1);

	/* average the (one-sided) power spectra */
	norm=1.0+decay*lastnorm;
	delta=decay*lastnorm/norm;
	for (i=0;i<npoint/2;i++) {
		im=(re=i+i)+1;
		spec=windata[re]*windata[re]+windata[im]*windata[im];

		if (initialize) {
			meanp[i]*=delta;
			meanp[i]+=spec/norm;
			deviation=spec-meanp[i];
			var=fabs(deviation);
			meanv[i]*=delta;
			meanv[i]+=var/norm;
		}
 		else {
			deviation=spec-meanp[i];
			var=fabs(deviation);

			freq=i*srate/npoint;
			/* if we are far from norm */
 
			if (var>THRESHOLD*meanv[i]) {
				fprintf(fpout,"%.2f\t%.2f\n",tstart,freq);
			}
			if (var<INCLUDE*meanv[i]) {
				meanp[i]*=delta;
				meanp[i]+=spec/norm;
				meanv[i]*=delta;
				meanv[i]+=var/norm;
			}
		}
	}

	if (initialize) initialize--;
	/* set averaging constant for next time */
	lastnorm=norm;

	fflush(fpout);
}
return 0;
}
