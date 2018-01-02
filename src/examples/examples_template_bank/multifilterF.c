/* GRASP: Copyright 1997,1998  Bruce Allen */
/* multifilter.c
This code is intended for machines where computation is cheap,
and communication is expensive.  The processsing is organized as
master/slaves (or manager/workers!).  The master process sends out data
chunks to individual slave processes.  These slave processes analyze
the data against all templates, then return the largest signal values
obtained for each template, along with other parameters like the time of
coalescense and the phase of coalescence.  They then get a new data chunk.
If STORE_TEMPLATES is set to 1, then the filters are computed once,
then stored internally by each slave.  This is the correct choice if each
slave has lots of fast memory available to it.  If STORE_TEMPLATES is set
to 0, then the slaves recompute the templates each time they use them.
This is the correct choice if each slave has only small amounts of fast
memory available.
*/


/* BROADCAST THE SAMPLE RATE TO THE DIFFERENT PROCESSES! */

#include "mpi.h"
#include "mpe.h"
#include "grasp.h"

#define NPOINT 65536       /* The size of our segments of data (6.5 secs)              */
#define FLO 120.0          /* The low frequency cutoff for filtering                   */
#define HSCALE 1.e21       /* A convenient scaling factor; results independent of it   */
#define MIN_INTO_LOCK 0.0  /* Number of minutes to skip into each locked section       */
#define CHIRPLEN 18000     /* length of longest allowed chirp                          */
#define MMIN 1.2           /* min mass object, solar masses                            */
#define MMAX 1.6           /* max mass object, solar masses                            */
#define DATA_SEGMENTS 60   /* largest number of data segments to process               */
#define NSIGNALS 11        /* number of signal values computed for each template       */
#define STORE_TEMPLATES 1  /* 0: slaves recompute templates. 1: slaves save templates. */

void shiftdata();          
void realft(float*,unsigned long,int);
int fill_buffer();

struct Saved {
   float tstart;
   int gauss;
};

short *datas;
struct fgetinput fgetinput;
struct fgetoutput fgetoutput;
int npoint,remain=0,needed,diff,gauss_test,num_sent=0;
float *twice_inv_noise,*htilde,*data,*mean_pow_spec,tstart;
float srate=9868.4208984375,decaytime,*response;
double norm,decay,datastart;
FILE *fpifo,*fpss,*fplock;

int main(int argc,char *argv[]) 
{

int *lchirppoints,num_stored;
float *ltc,*lch0tilde,*lch90tilde;
int myid,numprocs,i,j,maxi,impulseoff,*chirppoints,indices[8],num_templates;
int slave,more_data,temp_no,num_recv=0,namelen,completed=0,longest_template=0;
float *tc,m1,m2,*template_list,*sig_buffer,distance,snr_max,var,timeoff,timestart;
float n0,n90,inverse_distance_scale,*output90,*output0,*ch0tilde,*ch90tilde;
float lin0,lin90,varsplit,stats[8],gammq(float,float);
double prob;
FILE *fpout;
MPI_Status status;
char processor_name[MPI_MAX_PROCESSOR_NAME],logfile_name[64],name[64];
struct Scope Grid;
struct Saved *saveme;

/* start MPI, find number of processes, find process number */
MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
MPI_Comm_rank(MPI_COMM_WORLD,&myid);
MPI_Get_processor_name(processor_name,&namelen);
MPE_Init_log();

/* number of points to sample and fft (power of 2) */
needed=npoint=NPOINT;

/* Gravity wave signal (frequency domain) & twice inverse noise power */
htilde=(float *)malloc(sizeof(float)*npoint+sizeof(float)*(npoint/2+1));
twice_inv_noise=htilde+npoint;

/* Structure for saving information about data sent to slaves */
saveme=(struct Saved *)malloc(sizeof(struct Saved)*numprocs);

/* MASTER */
if (myid==0) {
   MPE_Describe_state(1,2,"Templates->Slaves","red:vlines3");
   MPE_Describe_state(3,4,"Data->Slaves","blue:gray3");
   MPE_Describe_state(5,6,"Master Receive","brown:light_gray");
   MPE_Describe_state(7,8,"Data->Master","yellow:dark_gray");
   MPE_Describe_state(9,10,"Slave Receive","orange:white");
   MPE_Describe_state(13,14,"Slaves<-templates","gray:black");
   MPE_Describe_state(15,16,"compute template","lavender:black");
   MPE_Describe_state(17,18,"real fft","lawn green:black");
   MPE_Describe_state(19,20,"correlate","purple:black");
   MPE_Describe_state(21,22,"orthonormalize","wheat:black");
   MPE_Describe_state(23,24,"likelyhood test","light sky blue:black");

   /* Set parameters for the inspiral search  */
   Grid.m_mn=MMIN;
   Grid.m_mx=MMAX;
   Grid.theta=0.964;
   Grid.dp=2*0.00213;
   Grid.dq=2*0.020;
   Grid.f_start=140.0;

   /* construct template set covering parameter space, m1 m2 storage */
   template_grid(&Grid);
   num_templates=Grid.n_tmplt;
   printf("The number of templates being used is %d\n",num_templates);
   template_list=(float *)malloc(sizeof(float)*2*num_templates);

   /* put mass values into an array */
    for (i=0;i<Grid.n_tmplt;i++) {
      template_list[2*i]=Grid.templates[i].m1;
      template_list[2*i+1]=Grid.templates[i].m2;
      printf("Mass values are m1 = %f  m2 = %f\n",Grid.templates[i].m1,Grid.templates[i].m2);
   }
   fflush(stdout);

   /* storage for returned signals (NSIGNALS per template) */
   sig_buffer=(float *)malloc(sizeof(float)*num_templates*NSIGNALS);

   /* broadcast templates */
   MPE_Log_event(1,myid,"send");
   MPI_Bcast(&num_templates,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(template_list,2*num_templates,MPI_FLOAT,0,MPI_COMM_WORLD);
   MPE_Log_event(2,myid,"sent");

   /* stores ADC data as short integers */
   datas=(short*)malloc(sizeof(short)*npoint);

   /* stores ADC data in time & freq domain, as floats */
   data=(float *)malloc(sizeof(float)*npoint);

   /* The response function (transfer function) of the interferometer */
   response=(float *)malloc(sizeof(float)*(npoint+2));

   /* The autoregressive-mean averaged noise power spectrum */
   mean_pow_spec=(float *)malloc(sizeof(float)*(npoint/2+1));

   /* Set up noise power spectrum and decay time */
   norm=0.0;
   clear(mean_pow_spec,npoint/2+1,1);
   decaytime=10.0*npoint/srate;
   decay=exp(-1.0*npoint/(srate*decaytime));

   /* stores ADC data as short integers */
   fgetinput.npoint=needed;
   fgetinput.nchan=1;
   fgetinput.files=framefiles;
   fgetinput.calibrate=1;
   fgetinput.chnames=(char **)malloc(fgetinput.nchan*sizeof(char *));
   fgetinput.locations=(short **)malloc(fgetinput.nchan*sizeof(short *));
   fgetoutput.npoint=(int *)malloc(fgetinput.nchan*sizeof(int));
   fgetoutput.ratios=(int *)malloc(fgetinput.nchan*sizeof(int));
   fgetinput.locations[0]=datas;

   /* channel name */
   fgetinput.chnames[0]="IFO_DMRO";

   /* set up channel names for different cases */
   if (NULL!=getenv("GRASP_REALTIME")) {
      /* don't care if locked */
      fgetinput.inlock=0;
   }
   else {
      /* only locked  */
      fgetinput.inlock=1;
   }

   /* while not finished, loop over slaves */
   for (slave=1;slave<numprocs;slave++) {
      if (fill_buffer()) {
         /* if new data exists, then send it (nonblocking?) */
         fprintf(stderr,"Master broadcasting data segment %d\n",num_sent+1);
         MPE_Log_event(3,myid,"send");
         MPI_Send(htilde,NPOINT+NPOINT/2+1,MPI_FLOAT,slave,++num_sent,MPI_COMM_WORLD);
         MPE_Log_event(4,myid,"sent");
         saveme[slave-1].gauss=gauss_test;
         saveme[slave-1].tstart=datastart;
         shiftdata();
      }
      else {
         /* tell remaining processes to exit */
         MPE_Log_event(3,myid,"send");
         MPI_Send(htilde,NPOINT+NPOINT/2+1,MPI_FLOAT,slave,0,MPI_COMM_WORLD);
         MPE_Log_event(4,myid,"sent");
      }
   }
      
   /* now loop, gathering answers, sending out more data */
   while (num_sent!=num_recv) {
      more_data=fill_buffer();

      /*  listen for answer */
      MPE_Log_event(5,myid,"receiving...");
      MPI_Recv(sig_buffer,NSIGNALS*num_templates,MPI_FLOAT,MPI_ANY_SOURCE,
         MPI_ANY_TAG,MPI_COMM_WORLD,&status);
      MPE_Log_event(6,myid,"received");
      num_recv++;

      /* store the answers... */
      sprintf(name,"signals.%05d",status.MPI_TAG-1);
      fpout=fopen(name,"w");
      if (fpout==NULL) {
         fprintf(stderr,"Multifilter: can't open output file %s\n",name);
         MPI_Finalize();
         return 1;
      }
      fprintf(fpout,"# Gaussian %d\n",saveme[status.MPI_SOURCE-1].gauss);
      fprintf(fpout,"# tstart %f\n",saveme[status.MPI_SOURCE-1].tstart);
      fprintf(fpout,"# snr   distance    phase0     phase90     maxi\
      impulseoff impulsetime  startinspiral coalesce    variance   prob\n");
      for (i=0;i<num_templates;i++) {
         for (j=0;j<NSIGNALS-1;j++)
            fprintf(fpout,"%g\t",sig_buffer[i*NSIGNALS+j]);
         fprintf(fpout,"%f\n",sig_buffer[i*NSIGNALS+j]);

         /* if data stream has no obvious outliers, and chirp prob is high, print */
         if (sig_buffer[i*NSIGNALS+10]>0.03 && saveme[status.MPI_SOURCE-1].gauss) {
            printf("POSSIBLE CHIRP: signal file %d, template %d, SNR = %f, prob = %f\n",
                   status.MPI_TAG-1,i,sig_buffer[i*NSIGNALS],sig_buffer[i*NSIGNALS+10]);
            fflush(stdout);
         }

      }
      fclose(fpout);

      /* if there is more data, send it off */
      if (more_data) {
         fprintf(stderr,"Master broadcasting data segment %d\n",num_sent+1);
         MPE_Log_event(3,myid,"send");
         MPI_Send(htilde,NPOINT+NPOINT/2+1,MPI_FLOAT,status.MPI_SOURCE,++num_sent,MPI_COMM_WORLD);
         MPE_Log_event(4,myid,"sent");
         saveme[status.MPI_SOURCE-1].gauss=gauss_test;
         saveme[status.MPI_SOURCE-1].tstart=datastart;
         shiftdata();
      }
      /* or else tell the process that it can pack up and go home */
      else {
         printf("Shutting down slave process %d\n",status.MPI_SOURCE);
         MPE_Log_event(3,myid,"send");
         MPI_Send(htilde,NPOINT+NPOINT/2+1,MPI_FLOAT,status.MPI_SOURCE,0,MPI_COMM_WORLD);
         MPE_Log_event(4,myid,"sent");
      }
   }

   /* when all the answers are in, print results */
   printf("This is the master - all answers are in!\n");
}

/* SLAVES */
else {
   printf("Slave %d (%s) just got started...\n",myid,processor_name);
   fflush(stdout);

   /* allocate storage space */
   /* Ouput of matched filters for phase0 and phase pi/2, in time domain, and temp storage */
   output0=(float *)malloc(sizeof(float)*npoint);
   output90=(float *)malloc(sizeof(float)*npoint);

   /* get the list of templates to use */
   MPE_Log_event(13,myid,"receiving...");
   MPI_Bcast(&num_templates,1,MPI_INT,0,MPI_COMM_WORLD);
   sig_buffer=(float *)malloc(sizeof(float)*num_templates*NSIGNALS);
   template_list=(float *)malloc(sizeof(float)*2*num_templates);
   MPI_Bcast(template_list,2*num_templates,MPI_FLOAT,0,MPI_COMM_WORLD);
   MPE_Log_event(14,myid,"received");
   printf("Slave %d (%s) just got template list...\n",myid,processor_name);
   fflush(stdout);

   /* Orthogonalized phase 0 and phase pi/2 chirps, in frequency domain */
   num_stored=STORE_TEMPLATES*(num_templates-1)+1;
   lch0tilde=(float *)malloc(sizeof(float)*npoint*num_stored);
   lch90tilde=(float *)malloc(sizeof(float)*npoint*num_stored);
   lchirppoints=(int *)malloc(sizeof(float)*num_stored);
   ltc=(float *)malloc(sizeof(float)*num_stored);

   if (lch0tilde==NULL || lch90tilde==NULL || lchirppoints==NULL || ltc==NULL) {
      fprintf(stderr,"Node %d on machine %s: could not malloc() memory!\n",
              myid,processor_name);
      MPI_Abort(MPI_COMM_WORLD,1);
   }

   /* now enter an infinite loop, waiting for new inputs */
   while (1) {
      /* listen for data, parameters from master */
      MPE_Log_event(9,myid,"receiving...");
      MPI_Recv(htilde,NPOINT+NPOINT/2+1,MPI_FLOAT,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
      MPE_Log_event(10,myid,"received");
      printf("Slave %d (%s) got htilde (and noise spectrum) for segment %d \n",
             myid,processor_name,status.MPI_TAG);
      fflush(stdout);

      /* if this is a termination message, we are done! */
      if (status.MPI_TAG==0) break;

      /* compute signals */
      for (temp_no=0;temp_no<num_templates;temp_no++) {

         ch0tilde=lch0tilde+npoint*temp_no*STORE_TEMPLATES;
         ch90tilde=lch90tilde+npoint*temp_no*STORE_TEMPLATES;
         chirppoints=lchirppoints+temp_no*STORE_TEMPLATES;
         tc=ltc+temp_no*STORE_TEMPLATES;

         /* Compute the template, and store it internally, if desired */
         if (completed!=num_templates) {
            /* manufacture two chirps (dimensionless strain at 1 Mpc distance) */
            m1=template_list[2*temp_no];
            m2=template_list[2*temp_no+1];

            MPE_Log_event(15,myid,"computing");
            make_filters(m1,m2,ch0tilde,ch90tilde,FLO,npoint,srate,chirppoints,tc,4000,4);
            MPE_Log_event(16,myid,"computed");

            if (*chirppoints>longest_template) longest_template=*chirppoints;

            if (*chirppoints>CHIRPLEN) {
               fprintf(stderr,"Chirp m1=%f m2=%f length %d too long!\n",m1,m2,
                  *chirppoints);
               fprintf(stderr,"Maximum allowed length is %d\n",CHIRPLEN);
               fprintf(stderr,"Please recompile with larger CHIRPLEN value\n");
               fflush(stderr);
               MPI_Abort(MPI_COMM_WORLD,1);
            }

            /* normalize the chirp template */
            /* normalization of next line comes from GRASP (5.6.3) and (5.6.4) */
            inverse_distance_scale=2.0*HSCALE*(TSOLAR*C_LIGHT/MPC);
            for (i=0;i<*chirppoints;i++){
               ch0tilde[i]*=inverse_distance_scale;
               ch90tilde[i]*=inverse_distance_scale;
            }

            /* and FFT the chirps */
            MPE_Log_event(17,myid,"starting fft");
            realft(ch0tilde-1,npoint,1);
            MPE_Log_event(18,myid,"ending fft");
            MPE_Log_event(17,myid,"starting fft");
            realft(ch90tilde-1,npoint,1);
            MPE_Log_event(18,myid,"ending fft");

            if (STORE_TEMPLATES) completed++;

            /* print out the length of the longest template */
            if (completed==num_templates)
               printf("Slave %d: templates completed.  Longest is %d points\n",
                      myid,longest_template);
            fflush(stdout);
         }  /* done computing the template */

         /* orthogonalize the chirps: we never modify ch0tilde, only ch90tilde */
         MPE_Log_event(21,myid,"starting");
         orthonormalize(ch0tilde,ch90tilde,twice_inv_noise,npoint,&n0,&n90);
         MPE_Log_event(22,myid,"done");

         /* distance scale Mpc for SNR=1 */
         distance=sqrt(1.0/(n0*n0)+1.0/(n90*n90));

         /* find the moment at which SNR is a maximum */   
         MPE_Log_event(19,myid,"searching");
         find_chirp(htilde,ch0tilde,ch90tilde,twice_inv_noise,n0,n90,output0,output90,
                       npoint,CHIRPLEN,&maxi,&snr_max,&lin0,&lin90,&var);
         MPE_Log_event(20,myid,"done");

         /* identify when an impulse would have caused observed filter output */
         impulseoff=(maxi+*chirppoints)%npoint;
         timeoff=impulseoff/srate;
         timestart=maxi/srate;

         /* collect interesting signals to return */
         sig_buffer[temp_no*NSIGNALS]=snr_max;
         sig_buffer[temp_no*NSIGNALS+1]=distance;
         sig_buffer[temp_no*NSIGNALS+2]=lin0;
         sig_buffer[temp_no*NSIGNALS+3]=lin90;
         sig_buffer[temp_no*NSIGNALS+4]=maxi;
         sig_buffer[temp_no*NSIGNALS+5]=impulseoff;
         sig_buffer[temp_no*NSIGNALS+6]=timeoff;
         sig_buffer[temp_no*NSIGNALS+7]=timestart;
         sig_buffer[temp_no*NSIGNALS+8]=timestart+*tc;
         sig_buffer[temp_no*NSIGNALS+9]=var;

         prob=0.0;
         if (snr_max>5.0) {
            MPE_Log_event(23,myid,"testing");
            varsplit=splitup_freq2(lin0*n0/sqrt(2.0),lin90*n90/sqrt(2.0),ch0tilde,
                             ch90tilde,2.0/(n0*n0),twice_inv_noise,npoint,maxi,8,
                             indices,stats,output0,htilde);
            prob=gammq(4.0,4.0*varsplit);
            MPE_Log_event(24,myid,"done");
         }
         sig_buffer[temp_no*NSIGNALS+10]=prob;

      }  /* end of loop over the templates */

      /* return signals to master */
      MPE_Log_event(7,myid,"send");
      MPI_Send(sig_buffer,NSIGNALS*num_templates,MPI_FLOAT,0,status.MPI_TAG,MPI_COMM_WORLD);
      MPE_Log_event(8,myid,"sent");

   } /* end of loop over the data */
}

/* both slaves and master exit here */
printf("%s preparing to shut down (process %d)\n",processor_name,myid);
sprintf(logfile_name,"multifilter.%d.%d.log",numprocs,DATA_SEGMENTS);
MPE_Finish_log(logfile_name);
MPI_Finalize();
printf("%s shutting down (process %d)\n",processor_name,myid);
return 0;
}

void shiftdata() {
   int i;

   /* shift ends of buffer to the start */
   needed=npoint-CHIRPLEN+1;
   for (i=0;i<CHIRPLEN-1;i++) datas[i]=datas[i+needed];

   return;
}


/* This routine gets the data set, overlapping the data buffer as needed */
int fill_buffer() {
   int i,code,reset=0;
   static double lastcalibtime=-1.0;
   static int num_sent=0;

   if (num_sent==DATA_SEGMENTS)
      return 0;

   topofloop:

   /* get number of points required */
   fgetinput.npoint=needed;
   fgetinput.locations[0]=datas+npoint-needed;
   fgetinput.seek=0;

   code=fget_ch(&fgetoutput,&fgetinput);
   datastart=fgetoutput.dt-(npoint-needed)/fgetoutput.srate;
   srate=fgetoutput.srate;

   /* if nothing left, return 0 */
   if (code==0) return 0;

   /* if new locked section, skip forward */
   if (code==1) {
      fgetinput.seek=1;
      fgetinput.npoint=fgetoutput.srate*MIN_INTO_LOCK*60.0-needed;
      code=fget_ch(&fgetoutput,&fgetinput);   
      if (code==0) return 0;

      /* number of points needed will be fulllength */
      needed=npoint;
      reset=1;
      goto topofloop;
   }

   /* if we have recalibrated, get the response function, and put in scaling factor */
   if (lastcalibtime<fgetoutput.tcalibrate) {
      lastcalibtime=fgetoutput.tcalibrate;
      GRnormalize(fgetoutput.fri,fgetoutput.frinum,npoint,srate,response);
      for (i=0;i<npoint+2;i++) response[i]*=HSCALE/ARMLENGTH_1994;
   }

   /* copy integer data into floats */
   for (i=0;i<npoint;i++) data[i]=datas[i];

   /* find the FFT of data*/
   realft(data-1,npoint,1);

   /* normalized delta-L/L tilde */
   product(htilde,data,response,npoint/2);

   /* update the inverse of the auto-regressive-mean power-spectrum */
   avg_inv_spec(FLO,srate,npoint,decay,&norm,htilde,mean_pow_spec,twice_inv_noise);

   /* see if the data has any obvious outliers */
   gauss_test=is_gaussian(datas,npoint,-2048,2047,0);

   num_sent++;
   if (reset==1)
      return 1;
   else
      return 2;
}
