/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"
#include "time.h"

/* number of graphs that are used */
int numtodo=1;

/* should we make a postscript file */
#define PRINT 1

/* should we make a jpeg file (needs postscript file option set) */
#define MAKEJPEG 1

/* should we display on the screen */
#define DISPLAY 0

/* number of sample points used */
#define NPOINT (256*64*32)

/* number of points to group in logarithmic bins */
#define DISPLAYSCALE 1.005

const float norm=((0.25/65536.0)*(1.0/400.0)*(1.0/(2.0*M_PI)));
double sqrt(double);
void realft(float*,unsigned long, int);
double remove_mean(float *,int);
void smooth(float *,int,int);

int main(int argc,char **argv) {
   void graphout(float,float,int);
   float tstart=1.e35,srate=1.e-30,tmin,tmax,dt;
   float *dataf[3],freq,norm2;
   double time=0.0;
   int i,j,seq=0,code,npoint=NPOINT;
   short *data[3];
   struct fgetinput fgetinput;
   struct fgetoutput fgetoutput;
   char filename[256],command1[256],command2[256];
   FILE *fp;
   time_t time1,time2;
   float f0,f1,avgfreq,ligostd;
   double total[3],totfreq;
   int bins;

   /* number of channels */
   fgetinput.nchan=3;

   /* source of files */
   fgetinput.files=framefiles;

   /* storage for channel names, data locations, points returned, ratios */
   fgetinput.chnames=(char **)malloc(fgetinput.nchan*sizeof(char *));
   fgetinput.locations=(short **)malloc(fgetinput.nchan*sizeof(short *));
   fgetoutput.npoint=(int *)malloc(fgetinput.nchan*sizeof(int));
   fgetoutput.ratios=(int *)malloc(fgetinput.nchan*sizeof(int));

   /* set up channel names, etc. for different cases */
   fgetinput.chnames[0]="H0:PEM-LVEA_SEISX";  /* seismometer X */
   fgetinput.chnames[1]="H0:PEM-LVEA_SEISY";  /* seismometer Y */
   fgetinput.chnames[2]="H0:PEM-LVEA_SEISZ";  /* seismometer Z */
   fgetinput.inlock=0;

   /* number of points to get */
   fgetinput.npoint=npoint;

   /* don't seek, we need the sample values! */
   fgetinput.seek=0;

   /* but we don't need calibration information */
   fgetinput.calibrate=0;

   /* allocate storage space for data */
   for (j=0;j<3;j++) {
      data[j]=(short *)malloc(sizeof(short)*npoint);
      dataf[j]=(float *)malloc(sizeof(float)*npoint);
      fgetinput.locations[j]=data[j];
   }

   /* now loop ... */
      while (numtodo-->0) {
         /* seek, or get the sample values */
         code=fget_ch(&fgetoutput,&fgetinput);

         /* elapsed time, sample rate */
         tstart=fgetoutput.dt;
         srate=fgetoutput.srate;
		
         /* if no data left, return */
         if (code==0) return 0;

	 /* copy data to floats */
	for (j=0;j<3;j++) {
		for (i=0;i<npoint;i++) {
			dataf[j][i]=data[j][i];
		}
		/* fft the data */
		remove_mean(dataf[j],npoint);
		realft(dataf[j]-1,npoint,1);


		/* normalize to get one-sided amplitude spectrum in meters/rHz */
		norm2=norm*sqrt(2.0/(npoint*srate));
		for (i=1;i<npoint/2;i++) {
			int ii,ir;
			ir=2*i;
			ii=ir+1;
			freq=i*fgetoutput.srate/npoint;
			dataf[j][i]=norm2*sqrt(dataf[j][ir]*dataf[j][ir]+dataf[j][ii]*dataf[j][ii])/freq;
		}
	}

	/* smooth */
	for (j=0;j<3;j++) {
		dataf[j][0]=0.0;
		smooth(dataf[j],npoint/2,8);
	}

        /* output */
	sprintf(filename,"Spec.%d",(int)fgetoutput.tstart_gps);

	/* make file containing data (and average over the bins) */
	fp=fopen(filename,"w");
	/* initialize variables */
	for (j=0;j<3;j++) total[j]=0.0;
	totfreq=f0=0.0;
	f1=0.01;
	bins=0;
	/* step through points in array */
	for (i=1;i<npoint/2;i++) {
		freq=i*fgetoutput.srate/npoint;
		/* if in the range of the current frequency bin, accumulate */
		if (f0<=freq && freq<f1) {
			for (j=0;j<3;j++) total[j]+=dataf[j][i];
			totfreq+=freq;
			bins++;
		}
		else {
			/* output if bin non-empty and freq range correct */
			if (bins && (avgfreq=totfreq/bins)>0.01 && avgfreq<100.0) {
				/* compute LIGO standard spectrum (m/rHz) */
				if (avgfreq<1.0) ligostd=1.e-9/(avgfreq*avgfreq*avgfreq);
				else if (avgfreq<10.0) ligostd=1.e-9;
				else ligostd=1.e-7/(avgfreq*avgfreq);
				fprintf(fp,"%f %e %e %e %e\n",avgfreq,total[0]/bins,total[1]/bins,total[2]/bins,ligostd);
			}
			/* set location of next bin, and accumulate current point */
			f0=freq;
			f1=f0*DISPLAYSCALE;
			bins=1;
			for (j=0;j<3;j++) total[j]=dataf[j][i];
			totfreq=freq;
		}
	}
	fclose(fp);

        time1=fgetoutput.tstart_gps;
        time2=time1-7*3600;

	/* make a file containing the right paramters */
	fp=fopen("seismic.params","w");
	fprintf(fp,"focus g0\n");
	fprintf(fp,"g0 type logxy\n");
	fprintf(fp,"autoscale\n");
	fprintf(fp,"title \"Ground motion (x: black y: red z: green LIGO:blue)\"\n");
/*	fprintf(fp,"subtitle \"UTC %.24s ",asctime(gpstime(&time1))); */
	time1=time2+npoint/srate;
	fprintf(fp,"subtitle \"Pacific Time %.19s to ",asctime(gpstime(&time2)));
	fprintf(fp,"%.24s\"\n",asctime(gpstime(&time1)));
	fprintf(fp,"xaxis label \"Frequency (Hz)\"\n");
	fprintf(fp,"yaxis label \"meters/sqrt(Hz)\"\n");
	fprintf(fp,"xaxis ticklabel prec 2\n");
	fprintf(fp,"xaxis  tick major grid on\n");
	fprintf(fp,"xaxis  tick minor grid on\n");
	fprintf(fp,"yaxis  tick major grid on\n");
	fprintf(fp,"yaxis  tick minor grid on\n");
	fprintf(fp,"xaxis  tick minor 1\n");
	fprintf(fp,"xaxis  tick major 1\n");
	fprintf(fp,"yaxis  tick minor 1\n");
	fprintf(fp,"yaxis  tick major 1\n");
	fprintf(fp,"xaxis  tick minor color 7\n");
	fprintf(fp,"xaxis  tick major color 7\n");
	fprintf(fp,"yaxis  tick minor color 7\n");
	fprintf(fp,"yaxis  tick major color 7\n");
	fprintf(fp,"world ymin 1.e-12\n");
	fprintf(fp,"world ymax 1.e-4\n");
	fclose(fp);

	printf("Created data file: %s\n",filename);
	if (DISPLAY) {
		sprintf(command1,"xmgr -nxy %s -param seismic.params &"
			,filename);
        	system(command1);
	}
	if (PRINT) {
		sprintf(command2,"grbatch -nxy %s -param seismic.params -printfile %s.ps -device 2"
		,filename,filename);
		system(command2);
		printf("Created postscript file: %s.ps\n",filename);

		if (MAKEJPEG) {
			sprintf(command2,"gs -sDEVICE=jpeg -r72 -q -dNOPAUSE -sOutputFile=%s.jpeg -- %s.ps"
			,filename,filename);
			system(command2);
			printf("Created jpeg file: %s.jpeg\n",filename);
                system("rm -f seismic.params");
		}
	}
   }
   return 0;
}

void smooth(float *data,int length,int avglen) {
        float norm;
        double total;
        int i,j,bot,top;
        float *work,*start;

        work=(float *)malloc(sizeof(float)*(length+avglen));
        start=work+avglen/2;
        for (i=0;i<avglen;i++)
                work[i]=data[0];
        for (i=length;i<length+avglen;i++)
                work[i]=data[length-1];
        memcpy(start,data,length*sizeof(float));
        for (i=0;i<length;i++) {
                total=0.0;
                bot=-avglen/2;
                top=avglen+bot;

                for (j=bot;j<top;j++)
                        total+=start[i+j];
                total/=(top-bot);
                data[i]=total;
        }
        free(work);
        return;
}
