/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"
#define NPOINT 4096

int main() {
   void graphout(int,float,float);
   float data[NPOINT],average[NPOINT/2],response[NPOINT+4];
   float spec,decaytime;
   float srate,tstart=0,freq,tlock=0.0;
   int i,j,code,npoint,ir,ii,reset=0,pass=0;
   short datas[NPOINT];
   double mod;
   struct fgetinput fgetinput;
   struct fgetoutput fgetoutput;

   /* number of channels needed is one */
   fgetinput.nchan=1;

   /* use utility function framefiles() to retrieve file names */
   fgetinput.files=framefiles;

   /* storage for channel names, data locations, points returned, ratios */
   fgetinput.chnames=(char **)malloc(fgetinput.nchan*sizeof(char *));
   fgetinput.locations=(short **)malloc(fgetinput.nchan*sizeof(short *));
   fgetoutput.npoint=(int *)malloc(fgetinput.nchan*sizeof(int));
   fgetoutput.ratios=(int *)malloc(fgetinput.nchan*sizeof(int));

   /* set up channel name */
   fgetinput.chnames[0]="IFO_DMRO";

   /* set up channel names for different cases */
   if (NULL!=getenv("GRASP_REALTIME")) {
      /* for Caltech 40-meter lab */
      fgetinput.inlock=0;
   }
   else {
      /* for Nov 1994 data set */
      fgetinput.inlock=1;
   }

   /* number of points to sample and fft (power of 2) */
   fgetinput.npoint=npoint=NPOINT;

   /* we do need the data, so don't seek */
   fgetinput.seek=0;

   /* do need calibration information */
   fgetinput.calibrate=1;

   /* where to put the data points */
   fgetinput.locations[0]=datas;

   /* set the decay time (sec) */
   decaytime=30.0;

   /* get data */
   while ((code=fget_ch(&fgetoutput,&fgetinput))) {
      tstart=fgetoutput.dt;
      srate=fgetoutput.srate;

       /* put data into floats */
      for (i=0;i<npoint;i++) data[i]=datas[i];

      /* use the swept-sine calibration (properly interpolated) to get R(f) */
      if (!pass++) GRnormalize(fgetoutput.fri,fgetoutput.frinum,npoint,srate,response);

      /* Reset if just locked */
      if (code==1) {
         reset=0;
         tlock=tstart;
         avg_spec(data,average,npoint,&reset,srate,decaytime,2,1);
      } else {

         /* track average power spectrum, with Welch windowing.  */
         avg_spec(data,average,npoint,&reset,srate,decaytime,2,1);

         /* loop over all frequencies except DC (j=0) & Nyquist (j=npoint/2) */
         for (j=1;j<npoint/2;j++) {
            /* subscripts of real, imaginary parts */
            ii=(ir=j+j)+1;
            /* frequency of the point */
            freq=srate*j/npoint;
            /* determine power spectrum in (meters/rHz) & print it */
            mod=response[ir]*response[ir]+response[ii]*response[ii];
            spec=sqrt(average[j]*mod);
            printf("%e\t%e\n",freq,spec);
         }
         /* print out useful things for xmgr program ... */
         graphout(0,tstart,tlock);
      }
   }
   return 0;
}


void graphout(int last,float time,float tlock) {
   static int count=0;
   printf("&\n");                            /* end of set marker             */
   /* first time we draw the plot */
   if (count++==0) {
      printf("@doublebuffer true\n");       /* keeps display from flashing    */
      printf("@focus off\n");               /* turn off the focus markers     */
      printf("@s0 color 2\n");              /* FFT is red                     */
      printf("@g0 type logxy\n");           /* set graph type to log-log      */
      printf("@autoscale \n");              /* autoscale FFT                  */
      printf("@world xmin %e\n",10.0);      /* set min x                      */
      printf("@world xmax %e\n",5000.0);    /* set max x                      */
      printf("@world ymin %e\n",1.e-19);    /* set min y                      */
      printf("@world ymax %e\n",1.e-9);     /* set max y                      */
      printf("@yaxis tick minor on\n");     /* turn on tick marks             */
      printf("@yaxis tick major on\n");     /* turn on tick marks             */
      printf("@yaxis tick minor 2\n");      /* turn on tick marks             */
      printf("@yaxis tick major 1\n");      /* turn on tick marks             */
      printf("@redraw \n");                 /* redraw graph                   */
      printf("@xaxis label \"f (Hz)\"\n");  /* FFT horizontal axis label      */
      printf("@yaxis label \"meters/rHz\"\n");  /* FFT vertical axis label    */
      printf("@title \"Calibrated IFO Spectrum\"\n");/* set title             */
      /* set subtitle     */
      printf("@subtitle \"%.2f sec since last lock. t = %.2f sec.\"\n",time-tlock,time);
      if (!last) printf("@kill s0\n");      /* kill graph; ready to read agai */
   }
   else {
      /* other times we redraw the plot */
      /* set subtitle     */
      printf("@subtitle \"%.2f sec since last lock. t = %.2f sec.\"\n",time-tlock,time);
      printf("@s0 color 2\n");              /* FFT is red                      */
      printf("@g0 type logxy\n");           /* set graph type to log-log       */
      printf("@world xmin %e\n",10.0);      /* set min x                       */
      printf("@world xmax %e\n",5000.0);    /* set max x                       */
      printf("@world ymin %e\n",1.e-19);    /* set min y                       */
      printf("@world ymax %e\n",1.e-9);     /* set max y                       */
      printf("@yaxis tick minor on\n");     /* turn on tick marks              */
      printf("@yaxis tick major on\n");     /* turn on tick marks              */
      printf("@yaxis tick minor 2\n");      /* turn on tick marks              */
      printf("@yaxis tick major 1\n");      /* turn on tick marks              */
      printf("@redraw\n");                  /* redraw the graph                */
      if (!last) printf("@kill s0\n");      /* kill graph, ready to read again */
   }
   return;
}
