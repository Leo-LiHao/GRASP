/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"

int main(int argc,char **argv) {
   void graphout(float,float,int);
   float tstart=1.e35,srate=1.e-30,tmin,tmax,dt;
   double time=0.0;
   int remain,i,seq=0,code,npoint=4096,seek;
   FILE *fp,*fplock;
   short *data;

      /* open the IFO output file and lock file in correct path */
      fp=grasp_open("GRASP_DATAPATH","channel.0","r");
      fplock=grasp_open("GRASP_DATAPATH","channel.10","r");

     /* allocate storage space for data */
     data=(short *)malloc(sizeof(short)*npoint);
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
               if (tstart<tmin-(npoint+20.)/srate) seek=1; else seek=0;
               code=get_data(fp,fplock,&tstart,npoint,data,&remain,&srate,seek);
               /* if no data left, return */
               if (code==0) return 0;
               /* we need to be outputting now... */
               if (tmin<=tstart){
                    for (i=0;i<npoint;i++) {
                         time=tstart+i/srate;
                         printf("%f\t%d\n",time,data[i]);
                    }
                    /* put out information for the graphing program */
                    graphout(tstart,tstart+npoint/srate,(argc==1 && time>=tmax));
               }
               /* if we are done with this interval, try next one */
               if (time>=tmax) break;
          }
     }
     /* close files and return */
     fclose(fp);
     return 0;
}
/* This routine is pipes output into the xmgr graphing program */
void graphout(float x1,float x2,int last) {
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
      printf("@subtitle \"IFO output %d\"\n",count++);/* set IFO subtitle        */
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
