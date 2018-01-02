/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"

int main(int argc,char **argv) {
   void graphout(float,float,int,char*);
   float tstart=1.e35,srate=1.e-30,tmin,tmax,dt;
   double time=0.0;
   int i,seq=0,code,npoint=4096;
   short *data;
   struct fgetinput fgetinput;
   struct fgetoutput fgetoutput;
   char *nstring;

   /* number of channels */
   fgetinput.nchan=1;

   /* number of samples per update */
   if (NULL!=(nstring=getenv("GRASP_NSAMPLE"))) {
      npoint=atoi(nstring);
      fprintf(stderr,"Display points set to environment variable: GRASP_NSAMPLE = %d\n",npoint);
   }

   /* source of files */
   fgetinput.files=framefiles;

   /* storage for channel names, data locations, points returned, ratios */
   fgetinput.chnames=(char **)malloc(fgetinput.nchan*sizeof(char *));
   fgetinput.locations=(short **)malloc(fgetinput.nchan*sizeof(short *));
   fgetoutput.npoint=(int *)malloc(fgetinput.nchan*sizeof(int));
   fgetoutput.ratios=(int *)malloc(fgetinput.nchan*sizeof(int));

   /* set up channel names, etc. for different cases */
   fgetinput.chnames[0]="IFO_DMRO";

   /* set up for different cases */
   if (NULL!=getenv("GRASP_CHANNEL")) {
      /* 40 meter lab */
      fgetinput.chnames[0]=getenv("GRASP_CHANNEL");
      fgetinput.inlock=0;
   }
   else {
      /* Nov 1994 data set */
      fgetinput.inlock=1;
   /* fgetinput.inlock=0; */
   }

   /* number of points to get */
   fgetinput.npoint=npoint;

   /* don't seek, we need the sample values! */
   fgetinput.seek=0;

   /* but we don't need calibration information */
   fgetinput.calibrate=0;

   /* allocate storage space for data */
   data=(short *)malloc(sizeof(short)*npoint);
   fgetinput.locations[0]=data;

   /* handle case where user has supplied t or dt arguments */
   if (argc==1) {
      tmin=-1.e30;
      dt=2.e30;
      argc=-1;
   }

   /* now loop ... */
   seq=argc;
   while (argc!=1) {
      /* get the next start time and dt */
      if (argc!=-1) {
         sscanf(argv[seq-argc+1],"%f",&tmin);
         sscanf(argv[seq-argc+2],"%f",&dt);
         argc-=2;
      }
      /* calculate the end of the observation interval, and get data */
      tmax=tmin+dt;
      while (1) {
         /* decide whether to skip (seek) or get sample values */
         if (tstart<tmin-(npoint+20.)/srate)
            fgetinput.seek=1;
         else
            fgetinput.seek=0;

         /* seek, or get the sample values */
         code=fget_ch(&fgetoutput,&fgetinput);

         /* elapsed time, sample rate */         tstart=fgetoutput.dt;
         srate=fgetoutput.srate;
		
         /* if no data left, return */
         if (code==0) return 0;

         /* we need to be outputting now... */
         if (tmin<=tstart){
            for (i=0;i<npoint;i++) {
               time=tstart+i/srate;
               printf("%f\t%d\n",time,data[i]);
            }

            /* put out information for the graphing program */
            graphout(tstart,tstart+npoint/srate,(argc==1 && time>=tmax),fgetinput.chnames[0]);
         }
         /* if we are done with this interval, try next one */
         if (time>=tmax) break;
      }
   }
   return 0;
}

/* This routine is pipes output into the xmgr graphing program */
void graphout(float x1,float x2,int last,char* ch_name) {
   static int count=0;
   printf("&\n");                            /* end of set marker             */
   /* first time we draw the plot */
   if (count==0) {
      printf("@doublebuffer true\n");       /* keeps display from flashing    */
      printf("@s0 color 3\n");              /* IFO graph is green             */
      printf("@view 0.1, 0.1, 0.9, 0.45\n"); /* set the viewport for IFO       */
      printf("@with g1\n");                 /* reset the current graph to FFT */
      printf("@view 0.1, 0.6, 0.9, 0.95\n");/* set the viewport FFT           */
      printf("@with g0\n");                 /* reset the current graph to IFO */
      printf("@world xmin %f\n",x1);        /* set min x                      */
      printf("@world xmax %f\n",x2);        /* set max x                      */
      printf("@autoscale\n");               /* autoscale first time through   */
      printf("@focus off\n");               /* turn off the focus markers     */
      printf("@xaxis label \"t (sec)\"\n"); /* IFO axis label                 */
      printf("@fft(s0, 1)\n");              /* compute the spectrum           */
      printf("@s1 color 2\n");              /* FFT is red                     */
      printf("@move g0.s1 to g1.s0\n");     /* move FFT to graph 1            */
      printf("@with g1\n");                 /* set the focus on FFT           */
      printf("@g1 type logy\n");            /* set FFT to log freq axis       */
      printf("@autoscale\n");               /* autoscale FFT                  */
      printf("@subtitle \"Spectrum\"\n");   /* set the subtitle               */
      printf("@xaxis label \"f (Hz)\"\n");  /* FFT axis label                 */
      printf("@with g0\n");                 /* reset the current graph IFO    */
      printf("@subtitle \"IFO output %d\"\n",count++);/* set IFO subtitle       */
      if (!last) printf("@kill s0\n");      /* kill IFO; ready to read again  */
   }
   else {
      /* other times we redraw the plot */
      printf("@s0 color 3\n");              /* set IFO green                   */
      printf("@fft(s0, 1)\n");              /* FFT it                          */
      printf("@s1 color 2\n");              /* set FFT red                     */
      printf("@move g0.s1 to g1.s0\n");     /* move FFT to graph 1             */
      printf("@subtitle \"%s %d\"\n",ch_name,count++);/* set IFO subtitle        */
      printf("@world xmin %f\n",x1);        /* set min x                       */
      printf("@world xmax %f\n",x2);        /* set max x                       */
      printf("@autoscale yaxes\n");         /* autoscale IFO                   */
      printf("@clear stack\n");             /* clear the stack                 */
      if (!last) printf("@kill s0\n");      /* kill IFO data                   */
      printf("@with g1\n");                 /* switch to FFT                   */
      printf("@g1 type logy\n");            /* set FFT to log freq axis       */
      printf("@clear stack\n");             /* clear stack                     */
      if (!last) printf("@kill s0\n");      /* kill FFT                        */
      printf("@with g0\n");                 /* ready to read IFO again         */
   }
   return;
}
