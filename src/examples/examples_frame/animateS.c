/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"

#define N_SLOW_CH 234
#define SLOW_CH_NAME_SIZE 40
short getslow(char *chanlist,char *channame,short *datalist);

int main(int argc,char **argv) {
   void graphout(float,float,int,char*);
   float tstart=1.e35,srate=1.e-30,tmin,tmax,dt;
   char slow_ch_name[SLOW_CH_NAME_SIZE];
   double time=0.0;
   int i,j,seq=0,code,nsec=64;
   short *data;
   struct fgetinput fgetinput;
   struct fgetoutput fgetoutput;
   float *slowdata;

   /* number of channels */
   fgetinput.nchan=1;
   if ( getenv("GRASP_CHANNEL") == NULL ) {
     fprintf(stderr,"Slow channel name must be specified by environment variable GRASP_CHANNEL\n");
     fprintf(stderr,"To obtain a list of names type animateS -l\n");
     exit(1);
   }
   strcpy(slow_ch_name,getenv("GRASP_CHANNEL"));


     /* source of files */
   fgetinput.files=framefiles;

   /* storage for channel names, data locations, points returned, ratios */
   fgetinput.chnames=(char **)malloc(fgetinput.nchan*sizeof(char *));
   fgetinput.locations=(short **)malloc(fgetinput.nchan*sizeof(short *));
   fgetoutput.npoint=(int *)malloc(fgetinput.nchan*sizeof(int));
   fgetoutput.ratios=(int *)malloc(fgetinput.nchan*sizeof(int));
   /* storage must be assigned for fgetinput.datatype as data is not in shorts */
   fgetinput.datatype=(char *)malloc(fgetinput.nchan*sizeof(char));
   fgetoutput.slow_names=(char *)malloc(N_SLOW_CH*SLOW_CH_NAME_SIZE*sizeof(char));
   slowdata=(float *)malloc(nsec*sizeof(float));


   /* set up channel names, etc. for different cases */
   fgetinput.chnames[0]="Slow";
     
   fgetinput.inlock=0;
 
   /* number of points to get */
   fgetinput.npoint=nsec*N_SLOW_CH;

   /* don't seek, we need the sample values! */
   fgetinput.seek=0;

   /* but we don't need calibration information */
   fgetinput.calibrate=0;

   /* allocate storage space for data */
   data=(short *)malloc(sizeof(float)*nsec*N_SLOW_CH);
   fgetinput.locations[0]=data;

 
   /* handle case where user just wants to see slow channel names */
   if (argc==2 && strcmp(argv[1],"-l") == 0 ) {
     fgetinput.seek=0;
     code=fget_ch(&fgetoutput,&fgetinput);
     for (i=0;i<N_SLOW_CH;i++) {
       printf("%s\n",fgetoutput.slow_names+SLOW_CH_NAME_SIZE*i);
     }
     exit(0);   
   }


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
         if (tstart<tmin-(nsec+20.)/srate)
            fgetinput.seek=1;
         else
            fgetinput.seek=0;

         /* seek, or get the sample values */
         code=fget_ch(&fgetoutput,&fgetinput);

	 for (i=0;i<N_SLOW_CH;i++) {
	   if ( strcmp(slow_ch_name,fgetoutput.slow_names+SLOW_CH_NAME_SIZE*i) == 0) break;
	 }
	 if (i == N_SLOW_CH) {
	   fprintf(stderr,"%s not found in Slow Channel\n",slow_ch_name);
	   fprintf(stderr,"To obtain a list of names type animateS -l\n");
	   exit(1);
	 }
  
	 for (j=0;j<nsec;j++) {
	   slowdata[j]=((float *)data)[i+N_SLOW_CH*j];
	 }

         /* elapsed time, sample rate */
         tstart=fgetoutput.dt;
         srate=fgetoutput.srate;
		
         /* if no data left, return */
         if (code==0) return 0;

         /* we need to be outputting now... */
         if (tmin<=tstart){
	     for (j=0;j<nsec;j++) {
               time=tstart+j;
               printf("%f\t%f\n",time,slowdata[j]);
	     }


	     /* put out information for the graphing program */
	     graphout(tstart,tstart+nsec/srate,(argc==1 && time>=tmax),slow_ch_name); 
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
      printf("@subtitle \"%s output %d\"\n",ch_name,count++);/* set Channel subtitle        */
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

#define N_SLOW_CH 234
#define SLOW_CH_NAME_SIZE 40
#include <string.h>
#include <stdio.h>

short getslow(char *chanlist,char *channame,short *datalist) {
	int i;
	for (i=0;i<N_SLOW_CH;i++) {
		if (strcmp(channame,chanlist+i*SLOW_CH_NAME_SIZE)==0)
			return datalist[i];
	}
	fprintf(stderr,"Did not find a channel named %s in frame.\n",channame);
	exit(1);
}
