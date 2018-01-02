/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"
/* macro to define the standard mathematical forms of mod */
#define MOD(X) ((X)>=0?((X)%num_display):(num_display-1+((X+1)%num_display)))  
#define FMOD2PI(X) ((X)>=0.0?(fmod((X),2.0*M_PI)):(2.0*M_PI+fmod((X),2.0*M_PI))) 

void medfit(float x[],float y[],int npoints,float *a ,float *b,float *dev);
void graphout(float,float,float,float);

int main() {
	short *datas;
	int i,num_points,num_win,num_freq,padded_length,max_lines,num_removed,remain,maxpass;
	int npass=1,top=0,estimate,nbins,padding_factor,num_display,nprint,index,new;
	int minbin,maxbin;
	float nwdt,*data,*mtap_spec_init,*mtap_spec_final,tstart,srate,creal,cimag;
	float *phase,phase1=0.0,amp1,phase2,*times,*linfitx,*linfity,binpreset;
	float displaytime,t1,*amp,dbin,ffit,intercept,slope,deviation,maxamp,displayamp=1.0;
	double time=0.0,fpreset,tinitial=0.0,offset;
	struct removed_lines *line_list;
	FILE *fpifo,*fplock;

	
	/* ------------------------ USER DEFINABLE ----------------------------*/
	/* data length, padded length (powers of 2!) */
	num_points=2048;
	padding_factor=8;

	/* your best guess for the line frequency you want to track */
	fpreset=582.395;

	/* set non-zero if you want us to estimate the best-fit frequency */
	estimate=0;

	/* number of (padded) frequency bins (either side) to search near fpreset */
	nbins=5;

	/* the number of phase/amplitudes to display */
	num_display=150;

	/* number of taper windows to use, and time-freq bandwidth */
	num_win=5;
	nwdt=3.0;

	/* the number of passes to make within the line removal algorithm */
	maxpass=1;

	/* num_points=2048; padding_factor=8;fpreset=582.395; */
	/* num_points=2048; padding_factor=4;fpreset=300.0; */
	num_points=2048; padding_factor=8;fpreset=582.395;

	/* --------------------- END OF USER DEFINABLE ----------------------------*/

	/* open the IFO output file and lock file */
	fpifo =grasp_open("GRASP_DATAPATH","channel.0","r");
	fplock=grasp_open("GRASP_DATAPATH","channel.10","r");

	padded_length=padding_factor*num_points;

	/* num frequencies including DC, Nyquist */
	num_freq=1+padded_length/2;

	/* max number of lines to report on */
	max_lines=64;

	/* allocate storage */
	datas=(short *)malloc(sizeof(short)*num_points);
	data=(float *)malloc(sizeof(float)*num_points);
	mtap_spec_init=(float *)malloc(sizeof(float)*num_freq);
	mtap_spec_final=(float *)malloc(sizeof(float)*num_freq);
	line_list=(struct removed_lines *)malloc(sizeof(struct removed_lines)*max_lines);
	amp=(float *)malloc(sizeof(float)*num_display);
	phase=(float *)malloc(sizeof(float)*num_display);
	times=(float *)malloc(sizeof(float)*num_display);
	linfitx=(float *)malloc(sizeof(float)*num_display);
	linfity=(float *)malloc(sizeof(float)*num_display);

	while (npass>0) {
		/* get a section of data... */
		new=get_data(fpifo,fplock,&tstart,num_points,datas,&remain,&srate,0);

		if (new==0)
			return 0;

		if (new==1) {
			fprintf(stderr,"\aTracker: New Locked Segment at time %f\n",tstart);
			ffit=fpreset;
			npass=1;
			tinitial=tstart;
			top=0;
			time=0.0;
		}
		else {
			/* Needed because timestamps only accurate to 1msec! */
			time+=(tstart-tinitial)/((double)npass-1.0);
		}

		/* the bin in which we expect to find frequency fpreset */
		binpreset=fpreset*2.0*num_freq/srate;
		minbin=binpreset-nbins;
		if (minbin<0) minbin=0;
		maxbin=binpreset+nbins;
		if (maxbin>num_freq) maxbin=num_freq;

		/* copy short data to float data */
		for (i=0;i<num_points;i++)
			data[i]=datas[i];

		/*  This line simulates the addition of an exponentially decaying sinusoid 
		for (i=0;i<num_points;i++)
			data[i]+=40.0*exp(-2.0*M_PI*310.0*(time+i/srate)/(2.0*100000.0))*
				sin(2.0*M_PI*310.0*(time+i/srate));
		*/

		/*  This line simulates the addition of an pure sinusoid 
		for (i=0;i<num_points;i++)
			data[i]+=200.0*sin(2.0*M_PI*310.0*(time+i/srate));
		*/

		/* remove the spectral lines from the data set */
		remove_spectral_lines(data,num_points,padded_length,nwdt,num_win,max_lines,maxpass,
			&num_removed,line_list,mtap_spec_init,mtap_spec_final,0,minbin,maxbin);

		/* if we fail to remove a line, amplitude set to zero, phase RETAINS PRIOR VALUE */
		amp1=0.0;  
		/* look in the list of removed lines for the right one */
		for (i=0;i<num_removed;i++) {

			/* the closest bin to our estimated frequency */
			dbin=binpreset-line_list[i].index;
			if (fabs(dbin)<=nbins) {
				/* use Taylor series to estime C at freq of interest */
 				creal=line_list[i].re+dbin*line_list[i].dcdbr+
					0.5*dbin*dbin*line_list[i].d2cdb2r;
				cimag=line_list[i].im+dbin*line_list[i].dcdbi+
					0.5*dbin*dbin*line_list[i].d2cdb2i;
				amp1=2.0*sqrt(creal*creal+cimag*cimag);
				phase1=atan2(cimag,creal)+2.0*M_PI*fmod(fpreset*time,1.0);
				break;
			}
		}

		/* save data in a circular buffers *[0..num_display-1]  */
		amp[top]=amp1;
		phase[top]=FMOD2PI(phase1);
		times[top]=time;

		/* how many values are we going to output to the graph? */
		nprint=(npass<num_display)?npass:num_display;

		/* cut out a piece for the linear fit */
		if (npass>=2) {

			/* adjust the phases to avoid boundary jumps */
			offset=0.0;
			index=MOD(top-nprint+1);
			linfitx[0]=times[index];
			linfity[0]=phase[index];
			for (i=1;i<nprint;i++) {
				index=MOD(top+i-nprint+1);
				linfitx[i]=times[index];	
				if (phase[index]-phase[MOD(index-1)]>M_PI)
					offset-=2.0*M_PI;
				else if (phase[index]-phase[MOD(index-1)]<-M_PI)
					offset+=2.0*M_PI;
				linfity[i]=phase[index]+offset;
			}

			/* do a robust linear fit */
			medfit(linfitx-1,linfity-1,nprint,&intercept,&slope,&deviation);

			/* now see what frequency the best fit corresponds to */
			ffit=fpreset-slope/(2.0*M_PI);

			/* if we are assuming a fixed frequency (not adapting) */
			if (!estimate) {
				slope=intercept=0.0;
			}

			/* print out amplitude if non-zero */
			maxamp=0.0;
			for (i=0;i<nprint;i++) {
			  	index=MOD(top+i-nprint+1);
				if (amp[index]>0.0)
					printf("%e\t%e\n",linfitx[i]+tinitial,amp[index]);
				else
					/* won't appear on the graph - out of bounds */
					printf("%e\t%f\n",linfitx[i]+tinitial,-1.0);
				if (amp[index]>maxamp) maxamp=amp[index];
			}
			/* separate data sets */
			printf("&\n");
			/* print out phase  if non-zero amplitude */
			for (i=0;i<nprint;i++) {
				phase2=linfity[i];
				phase2=FMOD2PI((phase2-slope*linfitx[i]-intercept));
				if (phase2>M_PI)
					phase2-=2.0*M_PI;
				phase2=(180.0/M_PI)*phase2;
			  	index=MOD(top+i-nprint+1);
				if (amp[index]>0.0)
					printf("%.8e\t%.8f\n",linfitx[i]+tinitial,phase2);
				else
					/* won't appear on the graph - out of bounds */
					printf("%.8e\t%f\n",linfitx[i]+tinitial,-500.0);
			}
		/* set up scale of the x-axis */
		t1=tinitial+linfitx[0];
		displaytime=num_display*(num_points/srate);
		/* set up scale of the amplitude graph y-axis */
		if (maxamp>0.9*displayamp) { 
			displayamp=1.3*maxamp;
			fprintf(stderr,"\aTracker: Line at %f Hz, amplitude just increased\n",fpreset);
		}
                else if (maxamp<0.4*displayamp && maxamp>0.0)
			displayamp=1.3*maxamp;		     

		graphout(t1,t1+displaytime,ffit,displayamp);
		fflush(stdout);
		}
	
			/* now display, then kill set */
          	npass++;
		top=MOD(top+1);
	}

	fclose(fpifo);
	fclose(fplock);

	return 0;
}


void graphout(float t1,float t2,float freq,float displayamp) {
	static int count=0;
	int xmaj,xmin;
	float ymaj=0.0,ymin=1.0;
	int amprec;

	xmin=(t2-t1)/10.0;
	xmaj=5*xmin;

	if (ymin<=displayamp/10.0)
		while (ymin<=displayamp/10.0) {
			ymin*=2.0;
			ymaj=4.0*ymin;
		}
	else
		while (ymin>displayamp/10.0) {
			ymin/=10.0;
			ymaj=5.0*ymin;
		}

	amprec=(int)log10(ymaj);
	if (amprec>1)
		amprec=0;
	else
	  	amprec=1-amprec;

	/* end of set marker */
	printf("&\n");

	if (count==0) {
		/* first time we draw the plot */
		printf("@doublebuffer true\n");
		printf("@focus off\n");
	}
	printf("@with g0\n");	
	printf("@move g0.s1 to g1.s0\n");
	printf("@title \"\\-Line Tracker\"\n");
	printf("@subtitle \"best estimate f=%f Hz\"\n",freq);
	printf("@s0 linestyle 0\n");
	printf("@s0 symbol color 4\n");
	printf("@s0 symbol 2\n");
	printf("@s0 symbol size 0.28\n");
	/*   printf("@s0 fill with color\n");  */
	printf("@s0 symbol fill 1\n");
	printf("@view 0.15, 0.53, 0.95, 0.90\n");
	/* set up x-axis for amplitude */ 
	printf("@world xmin %e\n",t1);
	printf("@world xmax %e\n",t2);
	printf("@xaxis tick major %d\n",xmaj);
	printf("@xaxis tick minor %d\n",xmin);
	printf("@xaxis ticklabel prec 1\n");
	printf("@xaxis ticklabel off\n");
	printf("@yaxis label \"\\-amplitude (ADC counts)\"\n");
	printf("@world ymin %e\n",0.0);
	printf("@world ymax %e\n",displayamp);
	printf("@yaxis tick major %e\n",ymaj);
	printf("@yaxis tick minor %e\n",ymin);
	if (amprec<4)
		printf("@yaxis ticklabel prec %d\n",amprec);
	else {
		printf("@yaxis ticklabel format general\n");
		printf("@yaxis ticklabel prec %d\n",1);
	}
	/* now do phase plot */
	printf("@with g1\n");
	printf("@s0 linestyle 0\n");
	printf("@s0 linewidth 0\n");
	printf("@s0 symbol color 2\n");
	printf("@s0 symbol 2\n");
	printf("@s0 symbol size 0.28\n");
	/*   printf("@s0 fill with color\n");  */
	printf("@s0 symbol fill 1\n");
	printf("@view 0.15, 0.1, 0.95, 0.47\n");
	/* set up x-axis for phase */ 
	printf("@world xmin %e\n",t1);
	printf("@world xmax %e\n",t2);
	printf("@xaxis tick major %d\n",xmaj);
	printf("@xaxis tick minor %d\n",xmin);
	printf("@xaxis ticklabel prec 1\n");
	printf("@xaxis label \"\\-time (sec)\"\n");
	/* set up y-axis for phase */
	printf("@world ymin %e\n",-180.0);
	printf("@world ymax %e\n",180.0);
	printf("@yaxis tick major 90\n");
	printf("@yaxis tick minor 45\n");
	printf("@yaxis ticklabel prec 0\n");
	printf("@yaxis label \"\\-phase (degrees)\"\n");
	printf("@xaxis label \"\\-time (sec)\"\n");
	/* draw plot */
	printf("@redraw\n");
	printf("@kill s0\n");
	printf("@with g0\n");
	printf("@kill s0\n");
      	count++;
	return;
}
