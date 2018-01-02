/* GRASP: Copyright 1997,1998  Bruce Allen */
/* binary_search.c
 * (based on multifilter.c)
 * JC: 25 July 1997, PRB: 19 August 1997, JC: 27 October 1997, 29 November 1997
 * BA: 11-13 April 1998
 */

#include <stdlib.h>
#include "binary_params.h"
#include "binary_string.h"
#include "mpi.h"
#include "grasp.h"
#include "assert.h"
#include <unistd.h>
#if(ISMPE)
#include "mpe.h"
#endif
static char *rcsid="$Id: binary_search.c,v 1.23 1999/07/07 18:04:43 ballen Exp $\n$Name: RELEASE_1_9_8 $";

#ifndef M_LN2
#define M_LN2       0.69314718055994530942      /* log e2 */
#endif

/* global variables for passing data from get_calibrated_data() */
double datastart;
float *n_inv_noise,*htilde,*pow_renorm,srate=9868.4208984375;
int num_templates,npoint=NPOINT,new_lock,gauss_test,insert_chirp=INSERT_CHIRP;

/* global variables for passing between master and slave in non mpi mode */
float *template_list,*sig_buffer;

/* MPI global variables */
int numprocs,myid,namelen;
char processor_name[MPI_MAX_PROCESSOR_NAME],logfile_name[256],*cmd_name;
void export_environment(int,const char*);
pid_t process_id;
#if (RENICE)
char commandstring[256];
int nicevalue;
#endif

/* define prototypes */
void realft(float *,unsigned long ,int);
void corr_coef(float *, float *, float *, int, float *, float *, float *);
/* the following routine is the Numerical Recipes routine select().
   However its name has been changed to NRselect.  See section 2.5.2 of the GRASP
   manual for details!
   */
float NRselect(unsigned long, unsigned long, float []);
int get_calibrated_data();

/* time-reversed filters */
void make_retlifs(float m1, float m2, float *ch1, float *ch2, float fstart,
                  int n, float srate, int *filled, float *t_coal,
                  int err_cd_sprs, int order);

/* information to save about segments of data */
struct Saved {
  double tstart;
  int gauss;
  int segmentno;
};

int main(int argc,char *argv[]) {
  /* prototypes */
  void master();
  void slave();
  
  /* Get the name of this program */
  cmd_name = argv[0];
  /* start MPI, find number of processes, find process number */
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  MPI_Get_processor_name(processor_name,&namelen);
  
  /* export any environment variables from master to slaves (repeat as many times as needed) */
  export_environment(myid,"GRASP_MFPATH");
  export_environment(myid,"GRASP_TEMPLATE");
  export_environment(myid,"GRASP_KILLSCRIPT");
#if (DBX)
  export_environment(myid,"DISPLAY");
#endif /* DBX */
  process_id=getpid();
#if (RENICE)
  /* set priority to 19 for slaves, 17 for master */
  if (myid)
    nicevalue=19;
  else
    nicevalue=17;
#if (!DBX)
  sprintf(commandstring,"renice %d -p %d",nicevalue,process_id);
  system(commandstring);
  printf("Processor %s (number %d) process number %d set NICE to %d\n",
	 processor_name,myid,process_id,nicevalue);
  fflush(stdout);
#endif /* !DBX */
#endif /* (RENICE) */
  
  /* create a file containing a kill script for this job */
  {
    int i;
    FILE *fpkill;
    char filemode[2];
    filemode[1]='\0';
    for (i=0;i<numprocs;i++) {
      MPI_Barrier(MPI_COMM_WORLD);
      if (i==myid) {
	filemode[0]=(myid==0)?'w':'a';
	fpkill=grasp_open("GRASP_KILLSCRIPT","killscript",filemode);
	if (myid==0) fprintf(fpkill,"#\n");
	fprintf(fpkill,"rsh %s kill -9 %d\n",processor_name,process_id);
    fflush(fpkill);
	fclose(fpkill);
      }
    }
  }
  
  /* if desired, pause to attach debugger  */
  /* if (myid==0) sleep(30); */
  
#if (DBX)
  sprintf(commandstring,"xterm -sb -e gdb %s %d &",cmd_name,process_id);
  if (myid==0) system(commandstring);
  sleep(15);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
#if(ISMPE)
  MPE_Init_log();
#endif /* (ISMPE) */
  
  /* Gravity wave signal (frequency domain),
     twice inverse noise power,
     and power renormalization factor */
  {
     int factor=1;
     factor=(myid==0)?1:NPERSLAVE;
     htilde = (float *)malloc(factor*(sizeof(float)*npoint
              +sizeof(float)*(npoint/2+1)+sizeof(float)));
  }
  n_inv_noise = htilde + npoint;
  pow_renorm = n_inv_noise + npoint/2 + 1;
  
  /* In the MPI version of the code, call the master or slave */
  if (myid==0)
     master();
  else
     slave();
  
  /* shutdown process */
  printf("%s preparing to shut down (process %d)\n",processor_name,myid);
  fflush(stdout);
#if(ISMPE)
  sprintf(logfile_name,"%s.%d.%d.log",cmd_name,numprocs,DATA_SEGMENTS);
  MPE_Finish_log(logfile_name);
#endif
  printf("%s waiting at MPI_Barrier... (process %d)\n",processor_name,myid);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  printf("%s shutting down (process %d)\n",processor_name,myid);
  fflush(stdout);
  return 0;
}

/* Function executed by the master node */
void master()
{
  MPI_Status status;
  struct Saved *saveme;
  int i,islave,num_sent=0,num_recv=0,*tot_per_slave,nperslave,*slave_shutdown;
  int startsegment=0;
  char *startsegenv;
  
  startsegenv=getenv("GRASP_STARTSEGMENT");
  if (startsegenv==NULL)
    startsegment=0;
  else {
    printf("Environment variable GRASP_STARTSEGMENT set to %s\n",startsegenv);
    startsegment=atoi(startsegenv);
    if (startsegment<0 || startsegment > 2048) startsegment=0;
    printf("Master starting with segment # %d\n",startsegment);
  }
  
#if(ISMPE)
  MPE_Describe_state(1,2,"Templates->Slaves","red:vlines3");
  MPE_Describe_state(3,4,"Data->Slaves","blue:gray3");
  MPE_Describe_state(5,6,"Master Receive","brown:light_gray");
  MPE_Describe_state(7,8,"Data->Master","yellow:dark_gray");
  MPE_Describe_state(9,10,"Slave Receive","orange:white");
  MPE_Describe_state(13,14,"Slaves<-templates","gray:black");
  MPE_Describe_state(15,16,"compute template","lavender:black");
  MPE_Describe_state(17,18,"real fft","lawn green:black");
  MPE_Describe_state(19,20,"correlate","purple:black");
  MPE_Describe_state(21,22,"correlation coefficients","wheat:black");
  MPE_Describe_state(23,24,"likelihood test","light sky blue:black");
#endif
  
  {  /* read in template list */
    
    float m1,m2;
    FILE *fptemplates;
    
    num_templates = NUM_TEMPLATES;
    template_list = (float *)malloc(sizeof(float)*2*num_templates);
    fptemplates=grasp_open("GRASP_TEMPLATE","templates.ascii","r");
    
    for (i=0;i<num_templates;i++) {
      if (EOF==fscanf(fptemplates,"%f %f\n",&m1,&m2)) {
        fprintf(stderr,"Warning: template file ended at template number %d\n%s\n",i,rcsid);
        fflush(stderr);
        num_templates = i;
        break;
      }
      template_list[2*i] = m1;
      template_list[2*i+1] = m2;
    }
    fclose(fptemplates);
  }
  
  /* print out the header */
  {
    time_t translate_time;
    float flo=FLO;
    int bytes=0,one=1;
    FILE *fpheader;
    
    fpheader = grasp_open("GRASP_MFPATH","signal.header","w");
    bytes += fprintf(fpheader,
		     "%d\n",HEADER_COMMENT_SIZE);
    bytes += fprintf(fpheader,
		     "This is output from a first pass filtering of the Nov 1994 data set.\n");
    time(&translate_time);
    bytes += fprintf(fpheader,"%s\n",ctime(&translate_time));
    bytes += fprintf(fpheader,"%s",comment);
    bytes += fprintf(fpheader,"%s",description);
    if (bytes>HEADER_COMMENT_SIZE) {
      fprintf(stderr,"HEADER_COMMENT_SIZE = %d must be increased to >= %d\n",HEADER_COMMENT_SIZE,bytes);
      fflush(stderr);
      MPI_Abort(MPI_COMM_WORLD,1);
    } 
    while (bytes<HEADER_COMMENT_SIZE) bytes += fprintf(fpheader,"\n");
    fwrite(&one,4,1,fpheader);
    fwrite(&num_templates,4,1,fpheader);
    fwrite(&flo,4,1,fpheader);
    fwrite(&srate,4,1,fpheader);
    for (i=0;i<2*num_templates;i++) fwrite(template_list+i,4,1,fpheader);
    fclose(fpheader);
  }
  
  /* storage for returned signals (NSIGNALS per template) */
  sig_buffer = (float *)malloc(sizeof(float)*num_templates*NSIGNALS*NPERSLAVE);
  
  /* Structure for saving information about data sent to slaves */
  saveme = (struct Saved *)malloc(sizeof(struct Saved)*numprocs*NPERSLAVE);
  tot_per_slave=(int *)malloc(sizeof(int)*numprocs);
  slave_shutdown=(int *)malloc(sizeof(int)*numprocs);
  { 
    int k;
    for (k=0;k<numprocs;k++) tot_per_slave[k]=slave_shutdown[k]=0;
  }
  
  /* broadcast templates */
#if(ISMPE)
  MPE_Log_event(1,myid,"send");
#endif
  MPI_Bcast(&num_templates,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(template_list,2*num_templates,MPI_FLOAT,0,MPI_COMM_WORLD);
#if(ISMPE)
  MPE_Log_event(2,myid,"sent");
#endif
  
  /* Skip segments 0 to startsegment-1 */
  {
    int i;
    i=startsegment;
    while(--i>=0) get_calibrated_data();
    num_sent=num_recv=startsegment;
  }
  
  /* while not finished, loop over slaves */
  for (islave=1;islave<numprocs;islave++)
    for (nperslave=0;(nperslave<NPERSLAVE) && (!slave_shutdown[islave]);nperslave++)
      if (get_calibrated_data()) { /* if new data exists, then send it (nonblocking?) */
	num_sent++;
	printf("Master broadcasting data segment %d to slave %d\n",num_sent,islave);
	fflush(stdout);
#if(ISMPE)
	MPE_Log_event(3,myid,"send");
#endif /* ISMPE */
	MPI_Send(htilde,NPOINT+NPOINT/2+1+1,MPI_FLOAT,islave,num_sent,MPI_COMM_WORLD);
#if(ISMPE)
	MPE_Log_event(4,myid,"sent");
#endif /* ISMPE */
	saveme[NPERSLAVE*(islave-1)+nperslave].gauss = gauss_test;
	saveme[NPERSLAVE*(islave-1)+nperslave].tstart = datastart;
	saveme[NPERSLAVE*(islave-1)+nperslave].segmentno = num_sent;
	tot_per_slave[islave]++;
      } else { /* tell remaining processes to exit */
	printf("Failed to get data segment %d\n",num_sent+1);
	fflush(stdout);
#if(ISMPE)
	MPE_Log_event(3,myid,"send");
#endif
	printf("Master - sent shutdown message to process %d\n",islave);
	fflush(stdout);
	MPI_Send(htilde,NPOINT+NPOINT/2+1+1,MPI_FLOAT,islave,0,MPI_COMM_WORLD);
	saveme[NPERSLAVE*(islave-1)+nperslave].segmentno=0;
	slave_shutdown[islave]=1;
#if(ISMPE)
	MPE_Log_event(4,myid,"sent");
#endif
      }
  
  /* now loop, gathering answers, sending out more data */
  printf("Entering loop to gather answers\n");
  fflush(stdout);
  while (num_sent!=num_recv) {
    
    FILE *fpout;
    char fname[256];
    int more_data;
    
    /*  listen for answer */
#if(ISMPE)
    MPE_Log_event(5,myid,"receiving...");
#endif
    
    /* This next block retrieves signals from the slaves */
    {
      int flag=0;
      while (flag==0) {
	/* check to see if a signal is waiting for pickup */
	MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&flag,&status);
	if (flag) {
	  /* if it is, then get it and exit */
	  MPI_Recv(sig_buffer,NSIGNALS*num_templates*NPERSLAVE,MPI_FLOAT,
		   MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
	}
	else
	  sleep(1);
      }
#if (DEBUG_COMM)
      printf("Master just got a response from slave %d (should be %d segments)\n",
	     status.MPI_SOURCE,tot_per_slave[status.MPI_SOURCE]);
      fflush(stdout);
#endif
      if (status.MPI_TAG != tot_per_slave[status.MPI_SOURCE]) {
	fprintf(stderr,"SERIOUS PROBLEM: slave node %d returned data for: %d data segments\n",
		status.MPI_SOURCE,status.MPI_TAG);
	fprintf(stderr,"BUT the master recorded that it was only sent:    %d data segments\n",
		tot_per_slave[status.MPI_SOURCE]);
	fflush(stderr);
      }
    }
#if(ISMPE)
    MPE_Log_event(6,myid,"received");
#endif /* ISMPE */
    num_recv+=tot_per_slave[status.MPI_SOURCE];
    
    /* Determine the correct file name */
    for (nperslave=0;nperslave<tot_per_slave[status.MPI_SOURCE];nperslave++) {
      
      sprintf(fname,"signal.%05d",
	      (saveme[NPERSLAVE*(status.MPI_SOURCE-1)+nperslave].segmentno)-1);
      
      /* open signal file for output */
      fpout=grasp_open("GRASP_MFPATH",fname,"w");
      
      /* Print the output data to a file */
      fwrite(&(saveme[NPERSLAVE*(status.MPI_SOURCE-1)+nperslave].tstart),8,1,fpout);
      fwrite(&(saveme[NPERSLAVE*(status.MPI_SOURCE-1)+nperslave].gauss), 4,1,fpout);
      fwrite(sig_buffer+nperslave*NSIGNALS*num_templates,NSIGNALS*4,num_templates,fpout);
      fclose(fpout);
    } /* end of loop over nperslave */
    tot_per_slave[status.MPI_SOURCE]=0;

    for (nperslave=0;nperslave<NPERSLAVE;nperslave++) {
      more_data = get_calibrated_data();
#if (DEBUG_COMM)
      printf("\t\t\t\t\t\t\tMORE DATA RETURNED %d\n",more_data);
      fflush(stdout);
#endif /* DEBUG COMM */
      if (more_data) { /* if there is more data, send it off */
	num_sent++;
	printf("Master broadcasting data segment %d to slave %d\n",num_sent,status.MPI_SOURCE);
	fflush(stdout);
#if(ISMPE)
	MPE_Log_event(3,myid,"send");
#endif /* ISMPE */
	MPI_Send(htilde,NPOINT+NPOINT/2+1+1,MPI_FLOAT,status.MPI_SOURCE,num_sent,MPI_COMM_WORLD);
#if(ISMPE)
	MPE_Log_event(4,myid,"sent");
#endif /* ISMPE */
	tot_per_slave[status.MPI_SOURCE]++;
	saveme[NPERSLAVE*(status.MPI_SOURCE-1)+nperslave].gauss = gauss_test;
	saveme[NPERSLAVE*(status.MPI_SOURCE-1)+nperslave].tstart = datastart;
	saveme[NPERSLAVE*(status.MPI_SOURCE-1)+nperslave].segmentno = num_sent;
      } else if (!slave_shutdown[status.MPI_SOURCE]) {
#if(ISMPE)
	MPE_Log_event(3,myid,"send");
#endif /* ISMPE */
	MPI_Send(htilde,NPOINT+NPOINT/2+1+1,MPI_FLOAT,status.MPI_SOURCE,0,MPI_COMM_WORLD);
	saveme[NPERSLAVE*(status.MPI_SOURCE-1)+nperslave].segmentno=0;
	slave_shutdown[status.MPI_SOURCE]=1;
	printf("Master - sent shutdown message to process %d\n",status.MPI_SOURCE);
	fflush(stdout);
#if(ISMPE)
	MPE_Log_event(4,myid,"sent");
#endif /* ISMPE */
      }
    }
    
  }  /* end of while (num_sent!=num_recv) loop */
  
  /* when all the answers are in, print results */
  printf("This is the master - all answers are in!\n");
  fflush(stdout);
  
  /* shut down any remaining slaves */
  for (islave=1;islave<numprocs;islave++)
    if (slave_shutdown[islave]==0) {
      MPI_Send(htilde,NPOINT+NPOINT/2+1+1,MPI_FLOAT,islave,0,MPI_COMM_WORLD);
      printf("Master - sent shutdown message to process %d (don't know why still running!)\n",
	     islave);
      fflush(stdout);
    }
  free(saveme);
  return;
}

void slave()
{
  void realft(float *, unsigned long, int);
  static float *output0=NULL,*output90=NULL,*snrdata=NULL;
  static float *ltc=NULL,*lch0tilde=NULL,*lch90tilde=NULL;
  static int *lchirppoints,completed=0;
  int i,num_stored,temp_no,nskip=4,validdata=0,nperslave,location;
  float *htildel,*n_inv_noisel,*pow_renorml;
  MPI_Status status;
  int segno[NPERSLAVE];
  
#if(PART_TEMPLATE>-1)
  char templatefile[256];
  FILE *fp_pt[NPERSLAVE];
#endif
  
  printf("Slave %d (%s) just got started...\n",myid,processor_name);
  fflush(stdout);
  
  /* allocate storage space */
  /* Ouput of matched filters for phase0 and phase pi/2, in time domain, and temp storage */
  snrdata = (float *)realloc(snrdata,sizeof(float)*npoint/nskip);
  output0 = (float *)realloc(output0,sizeof(float)*npoint);
  output90 = (float *)realloc(output90,sizeof(float)*npoint);
  
  /* get the list of templates to use */
#if(ISMPE)
  MPE_Log_event(13,myid,"receiving...");
#endif
  MPI_Bcast(&num_templates,1,MPI_INT,0,MPI_COMM_WORLD);
  sig_buffer = (float *)malloc(sizeof(float)*num_templates*NSIGNALS*NPERSLAVE);
  template_list = (float *)malloc(sizeof(float)*2*num_templates);
  MPI_Bcast(template_list,2*num_templates,MPI_FLOAT,0,MPI_COMM_WORLD);
#if(ISMPE)
  MPE_Log_event(14,myid,"received");
#endif
  printf("Slave %d (%s) just got template list...\n",myid,processor_name);
  fflush(stdout);
  
  /* Phase 0 and phase pi/2 chirps, in frequency domain */
  num_stored = STORE_TEMPLATES*(num_templates - 1) + 1;
  lch0tilde = (float *)realloc(lch0tilde,sizeof(float)*npoint*num_stored);
  lch90tilde = (float *)realloc(lch90tilde,sizeof(float)*npoint*num_stored);
  lchirppoints = (int *)realloc(lchirppoints,sizeof(float)*num_stored);
  ltc = (float *)realloc(ltc,sizeof(float)*num_stored);
  
  if (lch0tilde==NULL || lch90tilde==NULL || lchirppoints==NULL || ltc==NULL) {
    fprintf(stderr,"Node %d on machine %s: could not malloc() memory!\n%s\n",
	    myid,processor_name,rcsid);
    fflush(stderr);
    MPI_Abort(MPI_COMM_WORLD,5);
  }
  
  /* now enter an infinite loop, waiting for new inputs */
  while (1) {
    /* listen for data, parameters from master */
#if(ISMPE)
    MPE_Log_event(9,myid,"receiving...");
#endif
    validdata=0;
    for (nperslave=0;nperslave<NPERSLAVE;nperslave++) {
      MPI_Recv(htilde+nperslave*(NPOINT+NPOINT/2+1+1),NPOINT+NPOINT/2+1+1,MPI_FLOAT,0,MPI_ANY_TAG,
	       MPI_COMM_WORLD,&status);
      segno[nperslave]=status.MPI_TAG;
#if (DEBUG_COMM)
      printf("Slave %s node %d loop: message tag is %d\n",processor_name,myid,status.MPI_TAG);
      fflush(stdout);
#endif /* DEBUG COMM */
      if (segno[nperslave]!=0)
	validdata++;
      else
	break;
    }
#if(ISMPE)
    MPE_Log_event(10,myid,"received");
#endif
    /* if this is a termination message, we are done! */
    if (validdata==0)
      break;
    
    printf("Slave %d (%s) got htilde (and noise spectrum) for %d segments: %d to %d\n",
	   myid,processor_name,validdata,segno[0],segno[validdata-1]);
    fflush(stdout);
    
    /* compute signals */
    for (temp_no=0;temp_no<num_templates;temp_no++) {
      
      float *ch0tilde,*ch90tilde,*tc,r00,r11,r01,sigma2,distance;
      float snr_max,c0,c90,varsplit,stats[2*PBINS],med_renorm;
      int *chirppoints,num_crossing,maxi,impulseoff,indices[PBINS];
      
      /* To skip all templates except a particular one */
#if(PART_TEMPLATE>-1)
      for (nperslave=0;nperslave<validdata;nperslave++) {
	snr_max=0.0;
	varsplit=0.0;
	num_crossing=0;
	
	if((temp_no)!=PART_TEMPLATE) {
          sig_buffer[nperslave*num_templates*NSIGNALS+temp_no*NSIGNALS] = 0.0;
          sig_buffer[nperslave*num_templates*NSIGNALS+temp_no*NSIGNALS+1] = 0.0;
          sig_buffer[nperslave*num_templates*NSIGNALS+temp_no*NSIGNALS+2] = 0.0;
          sig_buffer[nperslave*num_templates*NSIGNALS+temp_no*NSIGNALS+3] = 0.0;
          ((int *)sig_buffer)[nperslave*num_templates*NSIGNALS+temp_no*NSIGNALS+4] = 0.0;
          sig_buffer[nperslave*num_templates*NSIGNALS+temp_no*NSIGNALS+5] = 0.0;          
          continue;  /* THIS IS NO LONGER RIGHT?? */
	}
	else{
          sprintf(templatefile,"template_%d.%d",PART_TEMPLATE,segno[nperslave]);
          fp_pt[nperslave] = grasp_open("GRASP_MFPATH",templatefile,"w");
	}
      }
#endif
      
      ch0tilde = lch0tilde + npoint*temp_no*STORE_TEMPLATES;
      ch90tilde = lch90tilde + npoint*temp_no*STORE_TEMPLATES;
      chirppoints = lchirppoints + temp_no*STORE_TEMPLATES;
      tc = ltc + temp_no*STORE_TEMPLATES;
      
      /* Compute the template, and store it internally, if desired */
      if (completed!=num_templates) {
        float m1,m2;
        int longest_template=0;
	
	/* manufacture two chirps (dimensionless strain at 1 Mpc distance) */
	m1 = template_list[2*temp_no];
	m2 = template_list[2*temp_no+1];
	
#if(ISMPE)
	if (!SMALL_MPELOG) MPE_Log_event(15,myid,"computing");
#endif
#if (REVERSE_FILTERS)
	make_retlifs(m1,m2,ch0tilde,ch90tilde,FLO,npoint,srate,chirppoints,tc,4000,4);
#else
	make_filters(m1,m2,ch0tilde,ch90tilde,FLO,npoint,srate,chirppoints,tc,4000,4);
#endif
	
#if(ISMPE)
	if (!SMALL_MPELOG) MPE_Log_event(16,myid,"computed");
#endif
	if (*chirppoints>longest_template) longest_template = *chirppoints;
	
	if (*chirppoints>CHIRPLEN) {
	  fprintf(stderr,"Chirp m1=%f m2=%f length %d too long!\n",m1,m2,
                  *chirppoints);
	  fprintf(stderr,"Maximum allowed length is %d\n",CHIRPLEN);
	  fprintf(stderr,"Please recompile with larger CHIRPLEN value\n");
          fprintf(stderr,"%s\n",rcsid);
	  fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD,6);      
	}
	
	/* normalize the chirp template */
	/* normalization of next line comes from GRASP (5.6.3) and (5.6.4) */
        {
	  float inverse_distance_scale=2.0*HSCALE*(TSOLAR*C_LIGHT/MPC);
	  for (i=0;i<*chirppoints;i++){
	    ch0tilde[i] *= inverse_distance_scale;
	    ch90tilde[i] *= inverse_distance_scale;
	  }
	}
	
	/* zero out the unused elements of the tilde arrays */
	for (i=(*chirppoints);i<npoint;i++)
	  ch0tilde[i]=ch90tilde[i]=0.0;
	
	/* and FFT the chirps */
#if(ISMPE)
	if (!SMALL_MPELOG) MPE_Log_event(17,myid,"starting fft");
#endif
	realft(ch0tilde-1,npoint,1);
#if(ISMPE)
	if (!SMALL_MPELOG) MPE_Log_event(18,myid,"ending fft");
	if (!SMALL_MPELOG) MPE_Log_event(17,myid,"starting fft");
#endif
	realft(ch90tilde-1,npoint,1);
#if(ISMPE)
	if (!SMALL_MPELOG) MPE_Log_event(18,myid,"ending fft");
#endif
	if (STORE_TEMPLATES) completed++;

	/* print out the length of the longest template */
	if (completed==num_templates)
	  fprintf(stderr,"Slave %d: templates completed.  Longest is %d points\n",
		  myid,longest_template);
	fflush(stdout);
      }  /* done computing the template */
      
      /* correlation coefficients */
      for (nperslave=0;nperslave<validdata;nperslave++) {
	n_inv_noisel=n_inv_noise+nperslave*(NPOINT+NPOINT/2+1+1);
	htildel=htilde+nperslave*(NPOINT+NPOINT/2+1+1);
	pow_renorml=pow_renorm+nperslave*(NPOINT+NPOINT/2+1+1);
#if(ISMPE)
	if (!SMALL_MPELOG) MPE_Log_event(21,myid,"starting");
#endif
	corr_coef(ch0tilde,ch90tilde,n_inv_noisel,npoint,&r00,&r11,&r01);
#if(ISMPE)
	if (!SMALL_MPELOG) MPE_Log_event(22,myid,"done");
#endif
	/* sigma squared and optimal distance scale Mpc for SNR=1            */
	sigma2 = 0.5*(r00 + r11);
	distance = sqrt(sigma2);
	
	/* find the moment at which SNR is a maximum */
#if(ISMPE)
	if (!SMALL_MPELOG) MPE_Log_event(19,myid,"searching");
#endif
	{
	  int count=0;
	  float x,y,amp,medsnr,expect=sqrt(M_LN2);
	  float x2py2,x2py2max,x2py2_thresh;
	  
	  correlate(output0,htildel,ch0tilde,n_inv_noisel,npoint);
	  correlate(output90,htildel,ch90tilde,n_inv_noisel,npoint);
	  
	  x2py2max=0.0;
	  num_crossing=0;
	  x2py2_thresh=THRESHOLD*THRESHOLD*sigma2;
	  maxi=-1;
	  for (i=PRESAFETY;i<NPOINT-POSTSAFETY;i++) {
	    x = output0[i];
	    y = output90[i];
            x2py2 = x*x + y*y;
#if(PART_TEMPLATE>-1)
	    fprintf(fp_pt[nperslave],"%d %f\n",i,sqrt(x2py2/sigma2));    
#endif
	    if (x2py2>x2py2max) {
	      x2py2max = x2py2;
	      maxi = i;
	    }
	    if (x2py2>x2py2_thresh) num_crossing++;
	    if (!(i%nskip)) snrdata[count++] = x2py2;
	  }
	  snr_max=sqrt(x2py2max/sigma2);
	  
#if(PART_TEMPLATE>-1)
	  fclose(fp_pt[nperslave]);
#endif
	  
	  /* check that we did correctly find a maximum offset */
	  assert(PRESAFETY <= maxi && maxi <= NPOINT-POSTSAFETY);
	  
	  c0 = (r11*output0[maxi] - r01*output90[maxi]);
	  c90 = (r00*output90[maxi] - r01*output0[maxi]);
	  amp = sqrt(c0*c0 + c90*c90);
	  
	  c0 /= amp;
	  c90 /= amp;
	  /* the following routine is the Numerical Recipes routine select().
	     However its name has been changed to NRselect.  See section 2.5.2 of the GRASP
	     manual for details!
	     */
	  medsnr = NRselect(count/2,count,snrdata-1);
	  med_renorm = expect/sqrt(medsnr/sigma2);
	}
	
#if(ISMPE)
	if (!SMALL_MPELOG) MPE_Log_event(20,myid,"done");
#endif
	/* identify when an impulse would have caused observed filter output */
	impulseoff = (maxi + *chirppoints)%npoint;
	
	/* collect interesting signals to return */
	location=NSIGNALS*num_templates*nperslave+temp_no*NSIGNALS;
	sig_buffer[location] = distance;
	sig_buffer[location+1] = snr_max;
	sig_buffer[location+2] = snr_max*(*pow_renorml);
	sig_buffer[location+3] = snr_max*med_renorm;
#if(COMPARE)
	((int *)sig_buffer)[location+4] = -maxi;
#else
	((int *)sig_buffer)[location+4] = impulseoff;     
#endif
	sig_buffer[location+5] = atan2(c90,c0);
	
	varsplit=0.0;
	if (snr_max>THRESHOLD) {
#if(ISMPE)
	  if (!SMALL_MPELOG) MPE_Log_event(23,myid,"testing");
#endif
#if (TWO_PHASE_R2)
          /* two-phase r^2 test */
	  varsplit = splitup_freq5(sqrt(0.5),sqrt(0.5),ch0tilde,ch90tilde,r00,
				   n_inv_noisel,npoint,
				   maxi,PBINS,indices,stats,output0,htildel);
#else
          /* Single phase r^2 test: */
	  varsplit = splitup_freq2(c0,c90,ch0tilde,ch90tilde,r00,
				   n_inv_noisel,npoint,
				   maxi,PBINS,indices,stats,output0,htildel);
#endif
	  varsplit /= sigma2;
#if(ISMPE)
	  if (!SMALL_MPELOG) MPE_Log_event(24,myid,"done");
#endif
	}
	sig_buffer[location+6] = varsplit;
	((int *)sig_buffer)[location+7] = num_crossing;
	
      }  /* end of loop over the templates */
    } /* end of loop over nperslave */
    
    /* return signals to master */
#if(ISMPE)
    MPE_Log_event(7,myid,"send");
#endif
#if (DEBUG_COMM)
    printf("Slave %s node %d sending message to master: tag is %d\n",processor_name,myid,
	   validdata);
    fflush(stdout);
#endif /* DEBUG_COMM */
    MPI_Send(sig_buffer,NSIGNALS*num_templates*NPERSLAVE,MPI_FLOAT,0,validdata,MPI_COMM_WORLD);
#if(ISMPE)
    MPE_Log_event(8,myid,"sent");
#endif
    if (validdata<NPERSLAVE) break;
  } /* end of loop over the data */
  
  printf("Node %d on machine %s: received a shutdown message\n",
	 myid,processor_name);
  fflush(stdout);
  return;
}




/* This routine exports environment variables to the nodes */
void export_environment(int node_number,const char *name) {
  char *environment,*var=NULL;
  int length;
  
  /* if this is the master, get enviroment variable */
  if (node_number==0) {
    if (name==NULL) {
      fprintf(stderr,"Argument to export_environment() can't be null!\n");
      fflush(stderr);
      MPI_Abort(MPI_COMM_WORLD,2);      
    }
    var=getenv(name);
    if (var==NULL) {
      fprintf(stderr,"Environment variable %s MUST be set!\n",name);
      fflush(stderr);
      MPI_Abort(MPI_COMM_WORLD,3);
    }
    /* calculate the length of the string to pass */
    length=strlen(name)+strlen(var)+2;
  }
  
  /* broadcast the length of the string */
  MPI_Bcast(&length,1,MPI_INT,0,MPI_COMM_WORLD);
  
  /* allocate storage for the string */
  if (NULL==(environment=(char *)malloc(length))) {
    fprintf(stderr,"Export_environment() node %d couldn't allocate memory!\n",
	    node_number);
    fflush(stderr);
    MPI_Abort(MPI_COMM_WORLD,4);
  }
  
  /*  print environment variable into string */
  if (node_number==0) {
    sprintf(environment,"%s=%s",name,var);
    fprintf(stderr,"Exporting environment variable %s\n",environment);
    fflush(stderr);
  }
  
  /* then broadcast it to the other processes */
  MPI_Bcast(environment,length,MPI_CHAR,0,MPI_COMM_WORLD);
  
  /* and put it into the local environment */
  if (node_number) {
    putenv(environment);
    fprintf(stderr,"Node %d receiving environment variable %s\n",
	    node_number,environment);
    fflush(stderr);
  }
  return;
}
