#include "grasp.h"
static char *rcsid="Development code";
#if defined (CLAPACK)
#include "f2c.h"
extern int chan_clean(int,int,float,float *,float *,float *,int,float **,float **,float *,float *,
		 complex *,complex *,float,complex *,complex *,integer,integer *);  
#else
typedef float real;                    /* from f2c.h */
typedef struct { real r, i; } complex; /* from f2c.h */
#endif

#define MIN(x,y) ((x) < (y) ? (x) : (y) ) 
#define MIN_BANDWIDTH 128
#define MAX_BANDWIDTH 128

void fopen_check(FILE *fp,char *fname);
char fft_dir[256],signal_name[256];

int main(int argc, char *argv[])
{
  int calc_rho(int,int,float,float *,float *,int,float **,float **,float *,
	       complex *,complex *,float *);
  void read_binary_fft(char *fname,int length,float *rp_fft,float *ip_fft,float *delta_f);
  void write_fft(char *fname,int length,float *rp_fft,float *ip_fft,float delta_f);
  void write_rho2(FILE *fp,float freq,float rho2,int nenv_chan,float *rho2_pairwise);
  void *errmalloc(char *arrayname,size_t bytes);
  void xmgr_files(int nenv_chan,char **chnames,int correlation_width);
  FILE *fp,*fp_rho2;
  char fname[256],rho2_fname[256],detector[256],cmd[256],**chnames,temp;
  complex *A=NULL, *B=NULL; 
  float **rp_env,**ip_env;
  float *rp_signal,*ip_signal,*rho2_pairwise,*rp_clean,*ip_clean;
  float threshold,rho2,modx2sum,delta_f,freq;
  int *chpts;
  int i,clean,offset,correlation_width,chan,n_chan,nenv_chan,length,check;
#if defined (CLAPACK)
  complex *R=NULL,*work=NULL;
  integer lwork,*ipivot=NULL;
#endif

  /*************************************************************************************/
  /*                       Read the configuration file                                 */
  /*************************************************************************************/
 
  if ( argc != 2) {
    printf("Usage: %s configuration-file-name\n",argv[0]);
    exit(1);
  }

  fp=fopen(argv[1],"r");
  fopen_check(fp,argv[1]);

  /*                       Comments start with a #   -                                  */
  /*    the first line after any comments should contain the detector name              */

  while (1) {
    fgets(detector,sizeof(detector),fp);
    if (detector[0] != '#') break;
  }  
  detector[strlen(detector)-1]='\0'; 

  /*    the next line gives the total number of channels including the `signal'         */ 
  check=fscanf(fp,"%d",&n_chan);
  if (check != 1) {
      fprintf(stderr,"Problems reading number of channels from %s\n",argv[1]); 
      exit(1);   
    }

  chnames=(char **)errmalloc("chnames",n_chan*sizeof(char *));
  for (i=0;i<n_chan;i++)   
     chnames[i]=(char *)malloc(256*sizeof(char));
  chpts=(int *)errmalloc("ip_signal",n_chan*sizeof(int));

  for (i=0;i<n_chan;i++) {
    check=fscanf(fp,"%s %c %d",chnames[i],&temp,&chpts[i]);
    /* the next n_chan lines of the configuration file should contain 3 columns         */
    /*                  - if not print an error message                                 */
    /*   the 3 columns are the channel name,the data type and the number of points      */
    /*     the data type is important for corr_init but is irrelevant here              */
    if (check != 3) {
      fprintf(stderr,"Problems reading 3-columns: channel-name data-type number-of-samples from %s\n",argv[1]); 
      exit(1);   
    }
  }

  /*          do we want to calculate (the fft of) the `cleaned' signal                 */ 
  /*    clean should be set to 1 if we want to calculate the cleaned signal             */
  check=fscanf(fp,"%d",&clean);
  if (check != 1) {
      fprintf(stderr,"Problems reading `clean' bit from %s\n",argv[1]); 
      exit(1);   
    }

  fclose(fp);

  /*************************************************************************************/
  /*           Determine the number of frequency bins we will use                      */ 
  /*************************************************************************************/

  length=chpts[0];
  for(i=1;i<n_chan;i++) 
     if (chpts[i]<length) length=chpts[i];

  length/=2;           /* Because length is the length of the complex FFT array */

  if ((length%MAX_BANDWIDTH) !=0) {
      length-=length%MAX_BANDWIDTH;
      fprintf(stderr,"Using %d FFT values (MAX_BANDWIDTH = %d)\n",length,MAX_BANDWIDTH);
  }

  strcpy(signal_name,chnames[0]);
  nenv_chan=n_chan-1;   /* because the configuration file includes the signal channel */ 

  /*************************************************************************************/
  /*                               Allocate memory                                     */		 
  /*************************************************************************************/

  /* Allocate memory for |x_i|^2,|y_i|^2 and x_i y_i*   */

  rp_signal=(float *)errmalloc("rp_signal",length*sizeof(float));
  ip_signal=(float *)errmalloc("ip_signal",length*sizeof(float));
  rp_env=(float **)errmalloc("rows of rp_env",nenv_chan*sizeof(float *));
  rp_env[0]=(float *)errmalloc("rows of rp_env[0]",nenv_chan*length*sizeof(float));
  for (i=1;i<nenv_chan;i++)
    rp_env[i]=rp_env[i-1]+length;

  ip_env=(float **)errmalloc("rows of ip_env",nenv_chan*sizeof(float *));
  ip_env[0]=(float *)errmalloc("rows of ip_env[0]",nenv_chan*length*sizeof(float));
  for (i=1;i<nenv_chan;i++)
    ip_env[i]=ip_env[i-1]+length;

  if (clean) { 
    rp_clean=(float *)errmalloc("rp_clean",length*sizeof(float));
    ip_clean=(float *)errmalloc("ip_clean",length*sizeof(float));
  }

  rho2_pairwise=(float *)errmalloc("rho2_pairwise",nenv_chan*sizeof(float));

  /*************************************************************************************/
  /*                      Allocate memory for lapack arrays                            */		 
  /*************************************************************************************/

  A=(complex *)errmalloc("A",nenv_chan*nenv_chan*sizeof(complex));
  B=(complex *)errmalloc("B",nenv_chan*sizeof(complex));
  if (clean) {
#if defined (CLAPACK)
    lwork=(integer) nenv_chan;
    R=(complex *)errmalloc("R",nenv_chan*sizeof(complex));
    work=(complex *)errmalloc("work",lwork*sizeof(complex));  
    ipivot=(integer *)errmalloc("ipivot",nenv_chan*sizeof(integer));
#endif
  } 

  /*************************************************************************************/
  /*************************************************************************************/
  /*                            Read data from file                                    */		 
  /*************************************************************************************/
  /*************************************************************************************/

  /*********************************************************************************/
  /*                       Determine name of data directory                        */		 
  /*********************************************************************************/
  for (i=0;i<256;i++) {
    if (argv[1][i] == '\0' || argv[1][i] == '.') {
      fft_dir[i] = '\0';
      break;
    }
    fft_dir[i]=argv[1][i];
  }
  strcat(fft_dir,"_fft");

  /*********************************************************************************/
  /*                       Read data from signal channel                           */		 
  /*********************************************************************************/
  sprintf(fname,"%s/%s-%s_fft.b",fft_dir,detector,chnames[0]);
  read_binary_fft(fname,length,rp_signal,ip_signal,&delta_f);

  /********************************************************************************/
  /*                 Read data from environmental channels                         */		 
  /********************************************************************************/
  for (chan=0;chan<nenv_chan;chan++) { 
    sprintf(fname,"%s/%s-%s_fft.b",fft_dir,detector,chnames[chan+1]);
    read_binary_fft(fname,length,rp_env[chan],ip_env[chan],&delta_f);
  } 

  /*************************************************************************************/
  /*                   Cycle through range of bandwidths                               */		 
  /*************************************************************************************/
  for (correlation_width=MIN_BANDWIDTH;correlation_width<=MAX_BANDWIDTH;correlation_width*=2) { 
  
    threshold=MIN(0.1,5.0/correlation_width); /* see discussion in Hua et al. */

    sprintf(rho2_fname,"%s/rho2_%s_%d.dat",fft_dir,signal_name,correlation_width);
    fp_rho2 = fopen(rho2_fname,"w");  
    fopen_check(fp_rho2,rho2_fname);
    fprintf(stderr,"Writing %s\n",rho2_fname); 

    /*************************************************************************************/
    /*                  Step through the range of offsets and call the                   */
    /*             function calc_rho where the major calculation is done                 */      	 
    /*************************************************************************************/

    for (offset=0;offset<length;offset+=correlation_width) { 

      calc_rho(offset,correlation_width,threshold,rp_signal,ip_signal,nenv_chan,
	       rp_env,ip_env,rho2_pairwise,A,B,&modx2sum);

     /*************************************************************************************/
     /*                     Calculate the `cleaned' signal                                */	
     /*************************************************************************************/

      rho2=0.0;                   /* set rho2=0.0 if we don't clean */
      if (clean) {
#if defined (CLAPACK)
	chan_clean(offset,correlation_width,threshold,&rho2,rp_signal,ip_signal,nenv_chan,
	 	   rp_env,ip_env,rp_clean,ip_clean,A,B,modx2sum,R,work,lwork,ipivot);	
#else
	fprintf(stderr,"Sorry cannot calculate cleaned channel without clapack installed\n"); 
	clean=0;
#endif
      } 


  /*************************************************************************************/
  /*    Print out |rho|^2 and the level of *pairwise* correlation between              */
  /*             the signal and each environmental channel in turn                     */       		 
  /*************************************************************************************/

     freq=(offset+0.5*correlation_width)*delta_f;
     write_rho2(fp_rho2,freq,rho2,nenv_chan,rho2_pairwise);
   
     
    }    /* end of offset loop */
 
    if (clean) {
      sprintf(fname,"%s/fftclean_%s_%d.dat",fft_dir,signal_name,correlation_width);
      write_fft(fname,length,rp_clean,ip_clean,delta_f);
    }

    /*************************************************************************************/
    /*                     Write xmgr parameter and shell files                          */		 
    /*************************************************************************************/

    xmgr_files(nenv_chan,chnames,correlation_width);
    sprintf(cmd,"chmod +x corr_view%d",correlation_width);
    system(cmd);
    sprintf(cmd,"corr_view%d &",correlation_width);
    system(cmd);


  } /* end of correlation_width loop */
 

  /*************************************************************************************/
  /*                               Free memory                                         */
  /*************************************************************************************/
  free(rp_signal);
  free(ip_signal);
  free(rp_env[0]);
  free(rp_env);
  free(ip_env[0]);
  free(ip_env);
  free(rp_clean);
  free(ip_clean);
  free(rho2_pairwise);
  /* free lapack arrays */
  free(A);
  free(B);
#if defined (CLAPACK)
  free(R);
  free(work);
  free(ipivot); 
#endif

  return(0);
}


void read_binary_fft(char *fname,int length,float *rp_fft,float *ip_fft,float *delta_f)
{
  /*************************************************************************************/
  /*    Read fft data from binary file: the first line contains the  frequency spacing */		 
  /*         then follow lines containing the real and imaginary part                  */
  /*************************************************************************************/
  FILE *fp;
  int i,check1,check2;

  fp = fopen(fname,"rb"); 
  fopen_check(fp,fname);
  fprintf(stderr,"Reading %s\n",fname); 
 
  /* first read  the frequency spacing */
  check1=fread(delta_f,sizeof(float),1,fp);
    if (check1 != 1) {
      fprintf(stderr,"Problems reading delta_f from %s\n",fname); 
      exit(1);   
    }  
  for (i=0;i<length;i++) {
    check1=fread((rp_fft+i),sizeof(float),1,fp);
    check2=fread((ip_fft+i),sizeof(float),1,fp);
   if ((check1 != 1) || (check2 != 1)) {
      fprintf(stderr,"Problems reading delta_f from %s\n",fname); 
      exit(1);   
    }  
  }

  fclose(fp); 
  rp_fft[0]=ip_fft[0]=0.0; /* set dc signal to 0 */  
}


void write_fft(char *fname,int length,float *rp_fft,float *ip_fft,float delta_f)
{
  /*************************************************************************************/
  /* Write fft data to an ascii file: the first line contains the  frequency spacing   */		 
  /*         then follow lines containing the real and imaginary part                  */
  /*************************************************************************************/  
  FILE *fp;
  int i;

  fp = fopen(fname,"w"); 
  fopen_check(fp,fname);
  fprintf(stderr,"Writing %s\n",fname);

  fprintf(fp,"%f\n",delta_f);
  /* We fill the DC component with data from the first frequency bin so that 
     we can do lin-log plots in xmgr without complaints.
     Note that the DC component is never used and is set to zero 
     by read_fft - still there should be a better way of doing this!  */

  fprintf(fp,"%f\t%f\n",rp_fft[1],ip_fft[1]);

  for (i=1;i<length;i++) {
    fprintf(fp,"%f\t%f\n",rp_fft[i],ip_fft[i]);
  }
  fclose(fp);
}


void write_rho2(FILE *fp,float freq,float rho2,int nenv_chan,float *rho2_pairwise)
{
  /*************************************************************************************/
  /*                     Write the list of rho values to file                          */		 
  /*************************************************************************************/
  int chan;

  fprintf(fp,"%.3f\t\t%.3f\t\t",freq,rho2); 
  for (chan=0;chan<nenv_chan;chan++) { 
    fprintf(fp,"%.3f ",rho2_pairwise[chan]); 
  }
  fprintf(fp,"\n"); 
}


void fopen_check(FILE *fp,char *fname) 
{
  /*************************************************************************************/
  /*            Checks to see if a file has been opened properly                       */
  /*             and if not write an appropriate error message                         */
  /*************************************************************************************/
  if ( fp == NULL ) {
    fprintf(stderr,"Problems opening %s\n",fname);
    exit(1);    
  }    
}


void *errmalloc(char *arrayname,size_t bytes) 
{
  /*************************************************************************************/
  /*            Allocate memory and print an error message if unsucessful              */
  /*************************************************************************************/ 
  void *pointer;
  pointer=malloc(bytes);
  if (pointer==NULL) {
    fprintf(stderr,"Cannot allocate %d bytes of memory for %s\n",(int) bytes,arrayname);
    exit(1);
  }
  return pointer;
}

void xmgr_files(int nenv_chan,char **chnames,int correlation_width)
{
  /*************************************************************************************/
  /*                     Write the xmgr parameter files                                */		 
  /*************************************************************************************/  
  FILE *fp;
  char param_fname[256],view_fname[256];
  int i;

  for (i=0;i<nenv_chan;i++) { 
  sprintf(param_fname,"xmgr.param%d_%d",correlation_width,i);

  fp = fopen(param_fname,"w"); 
  fopen_check(fp,param_fname);
  fprintf(stderr,"Writing %s\n",param_fname);

  fprintf(fp,"focus g%d\n",i);
  fprintf(fp,"autoscale\n");
  if (i == 0) {
    fprintf(fp,"xaxis label \"Correlations with %s: Frequency (Hz)\"\n",
	    chnames[0]);
    fprintf(fp,"xaxis  tick out\n");
    fprintf(fp,"xaxis  tick op bottom\n");
  }
  else {
    fprintf(fp,"xaxis  tick major off\n");
    fprintf(fp,"xaxis  tick minor off\n");
    fprintf(fp,"xaxis  ticklabel off\n");
    fprintf(fp,"yaxis  ticklabel start type spec\n");
    fprintf(fp,"yaxis  ticklabel start 0.5\n");
  }
  fprintf(fp,"yaxis  label layout perp\n");
  fprintf(fp,"yaxis  label char size 0.9\n");
  fprintf(fp,"yaxis  label place spec\n");
  fprintf(fp,"yaxis  label place -0.84, 0\n");
  fprintf(fp,"yaxis  ticklabel op right\n");

  fprintf(fp,"yaxis label \"%s\"\n",chnames[i+1]);  
  fprintf(fp,"yaxis  tick major 0.5\n");
  fprintf(fp,"yaxis  tick minor 0.25\n");
  fprintf(fp,"view ymin %f\n",0.1+0.88*i/nenv_chan);
  fprintf(fp,"view ymax %f\n",0.1+0.88*(i+1)/nenv_chan);    
  fprintf(fp,"world ymin 0\n");
  fprintf(fp,"world ymax 1\n"); 

  fprintf(fp,"s0 color %d\n",i+2);
  fprintf(fp,"\n");
  fclose(fp);
 }


  /*************************************************************************************/
  /*                     Write the xmgr shell files                                    */		 
  /*************************************************************************************/

  sprintf(view_fname,"corr_view%d",correlation_width);
  fp = fopen(view_fname,"w"); 
  fopen_check(fp,view_fname);
  fprintf(stderr,"Writing %s\n",view_fname);

  fprintf(fp,"xmgr -block %s/rho2_%s_%d.dat ",
	        fft_dir,signal_name,correlation_width);
    
  for (i=0;i<nenv_chan;i++) { 
    sprintf(param_fname,"xmgr.param%d_%d",correlation_width,i);
    fprintf(fp,"-graph %d -bxy 1:%d -param %s ",i,i+3,param_fname);
  }
  fprintf(fp,"&\n");
    
  fclose(fp);
}
