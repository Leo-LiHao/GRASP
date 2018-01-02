/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"
#define NPOINT 2048        /* determines frequency resolution */
#define MIN_INTO_LOCK 3.0  /* how many minutes to skip into locked sections */
#define SHIFT 1024         /* determines time resolution */
#define CHANNEL 0          /* channel number */        
#define AVG_TIME 10.0      /* seconds */
#define THRESHOLD 6.0      /* for outputting a point */
#define INCLUDE 10.0       /* threshold for including in averaging */


int main(int argc,char **argv) {
	void realft(float*,unsigned long,int);
	void clear(float*,int,int);
	FILE *fpifo,*fplock,*fpout=NULL;
	short *datas;
	int i,code,npoint,remain,im,re,nskip,pointsneeded,diff,initialize=10,filenumber=0;
	float spec,srate,tstart,*meanp,*meanv,*data,decaytime,datastart;
	float win,freq,deviation,var,*window;
	double decay=0.0,delta,norm,lastnorm=0.0;
	char name[64];

	/* open the IFO output file and lock file */
	sprintf(name,"channel.%d",CHANNEL);
	fpifo=grasp_open("GRASP_DATAPATH",name,"r");
	fplock=grasp_open("GRASP_DATAPATH","channel.10","r");

	/* number of points to sample and fft (power of 2) */
	npoint=pointsneeded=NPOINT;
	diff=npoint-pointsneeded;

	/* create arrays */
	datas=(short*)malloc(sizeof(short)*npoint);
	data=(float *)malloc(sizeof(float)*npoint);
	window=(float *)malloc(sizeof(float)*npoint);
	meanp=(float *)malloc(sizeof(float)*npoint/2);
	meanv=(float *)malloc(sizeof(float)*npoint/2);

	/* compute Welch window, once */
	for (i=0;i<npoint;i++) {
		win=(2.0*i)/npoint-1.0;
		window[i]=1.0-win*win;
	}

/* This is the main program loop, which aquires data, then filters it */
while ((code=get_data(fpifo,fplock,&tstart,pointsneeded,datas+diff,&remain,&srate,0))) {
	datastart=tstart-diff/srate;

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

		if (fpout!=NULL) fclose(fpout);
		sprintf(name,"ch%ddiag.%03d",CHANNEL,filenumber);
		fpout=fopen(name,"w");
 		fprintf(fpout,"@ sets linestyle 0\n");
 		fprintf(fpout,"@ sets symbol 1\n");
 		fprintf(fpout,"@ s%d symbol color %d\n",CHANNEL,CHANNEL+1);
		fprintf(fpout,"@ title \"19 November 1994 run 1\"\n");
		fprintf(fpout,"@ subtitle \"Time/Frequency statistics for channel %d\"\n",CHANNEL);
		fprintf(fpout,"@ xaxis label \"Time (sec)\"\n");
		fprintf(fpout,"@ yaxis label \"Frequency (Hz)\"\n");
		filenumber++;
 	}

	/* if just entering a locked stretch */
	if (code==1 && (nskip=MIN_INTO_LOCK*60*srate-pointsneeded)>0) {

		/* get three minutes total (seek)... */
		get_data(fpifo,fplock,&tstart,nskip,datas+diff,&remain,&srate,1);
		datastart=tstart-diff/srate;
 	}
	else {
		/* we are MIN_INTO_LOCK minutes into a locked stretch of data */		
		/* copy data into floats with windowing */
		for (i=0;i<npoint;i++)
			data[i]=window[i]*datas[i];

		/* find the FFT */
		realft(data-1,npoint,1);

		/* average the (one-sided) power spectra */
		norm=1.0+decay*lastnorm;
		delta=decay*lastnorm/norm;
		for (i=0;i<npoint/2;i++) {
			im=(re=i+i)+1;
			spec=data[re]*data[re]+data[im]*data[im];

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
					fprintf(fpout,"%.2f\t%.2f\n",datastart,freq);
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

		/* decide how many points are needed... */
		pointsneeded=SHIFT;
		if (remain<pointsneeded)
			pointsneeded=npoint;
		else {
			/* shift ends of buffer to the start */
			for (i=SHIFT;i<npoint;i++)
				datas[i-SHIFT]=datas[i];
			/* shift time marker for datas[] array (un-necessary?) */
			datastart+=SHIFT/srate;
		}
		diff=npoint-pointsneeded;
	}
}
return 0;
}
