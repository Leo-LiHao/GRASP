/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"
#define NPOINT 131072      /* The size of our segments of data (13.1 secs) */
#define FLO 120.0          /* The low frequency cutoff for filtering */
#define HSCALE 1.e21       /* A convenient scaling factor; results independent of it */
#define MIN_INTO_LOCK 3.0  /* Number of minutes to skip into each locked section */
#define SAFETY 1000        /* Padding safety factor to avoid wraparound errors */
#define DATA_SEGMENTS 3000 /* maximum number of data segments to filter */

struct fgetinput fgetinput;
struct fgetoutput fgetoutput;
short *datas;
int needed,npoint;
double datastart;
float srate=9868.4208984375,*response,min_into_lock=MIN_INTO_LOCK;

int main() {
   void realft(float*,unsigned long,int);
   int i,code,maxi,chirplen,impulseoff,chirppoints,indices[8],fill_buffer();
   float distance,snr_max,var,*mean_pow_spec;
   float *data,*htilde,*output90,*output0,*chirp0,*chirp90,*ch0tilde,*ch90tilde;
   float n0,n90,inverse_distance_scale,decaytime,*twice_inv_noise,tc;
   float lin0,lin90,invMpc_inject,varsplit,stats[8],gammq(float,float);
   double decay=0.0,norm,prob,timeoff,timestart;

   /* number of points to sample and fft (power of 2) */
   fgetinput.npoint=needed=npoint=NPOINT;

   /* stores ADC data as short integers */
   fgetinput.nchan=1;
   fgetinput.files=framefiles;
   fgetinput.calibrate=1;
   fgetinput.chnames=(char **)malloc(fgetinput.nchan*sizeof(char *));
   fgetinput.locations=(short **)malloc(fgetinput.nchan*sizeof(short *));
   fgetoutput.npoint=(int *)malloc(fgetinput.nchan*sizeof(int));
   fgetoutput.ratios=(int *)malloc(fgetinput.nchan*sizeof(int));
   fgetoutput.lastlock=fgetoutput.tstart=0.0;

   fgetinput.locations[0]=datas=(short*)malloc(sizeof(short)*npoint);

   /* set up channel names */
   fgetinput.chnames[0]="IFO_DMRO";

   /* set up channel names for different cases */
   if (NULL!=getenv("GRASP_REALTIME")) {
      /* don't care if locked */
      fgetinput.inlock=0;
      min_into_lock=0.0;
   }
   else {
      /* only locked  */
      fgetinput.inlock=1;
   }

   /* stores ADC data in time & freq domain, as floats */
   data=(float *)malloc(sizeof(float)*npoint);

   /* The phase 0 and phase pi/2 chirps, in time domain */
   chirp0=(float *)malloc(sizeof(float)*npoint);
   chirp90=(float *)malloc(sizeof(float)*npoint);

   /* Orthogonalized phase 0 and phase pi/2 chirps, in frequency domain */
   ch0tilde=(float *)malloc(sizeof(float)*npoint);
   ch90tilde=(float *)malloc(sizeof(float)*npoint);

   /* The response function (transfer function) of the interferometer */
   response=(float *)malloc(sizeof(float)*(npoint+2));

   /* The gravity wave signal, in the frequency domain */
   htilde=(float *)malloc(sizeof(float)*npoint);

   /* The autoregressive-mean averaged noise power spectrum */
   mean_pow_spec=(float *)malloc(sizeof(float)*(npoint/2+1));

   /* Twice the inverse of the mean noise power spectrum */
   twice_inv_noise=(float *)malloc(sizeof(float)*(npoint/2+1));

   /* Ouput of matched filters for phase0 and phase pi/2, in time domain, and temp storage */
   output0=(float *)malloc(sizeof(float)*npoint);
   output90=(float *)malloc(sizeof(float)*npoint);

   /* manufacture two chirps (dimensionless strain at 1 Mpc distance) */
   make_filters(1.4,1.4,chirp0,chirp90,FLO,npoint,srate,&chirppoints,&tc,0,4);
   /* normalization of next line comes from GRASP (5.6.3) and (5.6.4) */
   inverse_distance_scale=2.0*HSCALE*(TSOLAR*C_LIGHT/MPC);
   for (i=0;i<chirppoints;i++){
      ch0tilde[i]=chirp0[i]*=inverse_distance_scale;
      ch90tilde[i]=chirp90[i]*=inverse_distance_scale;
   }

   /* zero out the unused elements of the tilde arrays */
   for (i=chirppoints;i<npoint;i++)
      ch0tilde[i]=ch90tilde[i]=0.0;

   /* and FFT the chirps */
   realft(ch0tilde-1,npoint,1);
   realft(ch90tilde-1,npoint,1);

   /* set length of template including a safety margin */
   chirplen=chirppoints+SAFETY;
   if (chirplen>npoint) abort();

   /* This is the main program loop, which aquires data, then filters it */
   while (1) {
      /* fill the buffer with the number of points needed */
      code=fill_buffer();

      /* if no points left, we are done! */
      if (code==0) break;

      /* if just entering a new locked stretch, reset averaging over power spectrum */
      if (code==1) {
         norm=0.0;
         clear(mean_pow_spec,npoint/2+1,1);

         /* decay time for spectrum, in sec.  Set to 15x length of npoint sample */
         decaytime=15.0*npoint/srate;
         decay=exp(-1.0*npoint/(srate*decaytime));
      }

      /* copy integer data into floats */
      for (i=0;i<npoint;i++) data[i]=datas[i];

      /* inject signal in time domain (note output0[] used as temp storage only) */
      invMpc_inject=0.0;   /* To inject a signal at 10 kpc, set this to 100.0 */
      time_inject_chirp(1.0,0.0,12345,invMpc_inject,chirp0,chirp90,data,response,output0,npoint);

      /* find the FFT of data*/
      realft(data-1,npoint,1);

      /* normalized delta-L/L tilde */
      product(htilde,data,response,npoint/2);

      /* inject a signal in frequency domain, if desired */
      invMpc_inject=0.0;   /* To inject a signal at 10 kpc, set this to 100.0 */
      freq_inject_chirp(-0.406,0.9135,23456,invMpc_inject,ch0tilde,ch90tilde,htilde,npoint);

      /* update the inverse of the auto-regressive-mean power-spectrum */
      avg_inv_spec(FLO,srate,npoint,decay,&norm,htilde,mean_pow_spec,twice_inv_noise);

      /* orthogonalize the chirps: we never modify ch0tilde, only ch90tilde */
      orthonormalize(ch0tilde,ch90tilde,twice_inv_noise,npoint,&n0,&n90);
  
      /* distance scale Mpc for SNR=1 */
      distance=sqrt(1.0/(n0*n0)+1.0/(n90*n90));

      /* find the moment at which SNR is a maximum */
      find_chirp(htilde,ch0tilde,ch90tilde,twice_inv_noise,n0,n90,output0,output90,
                 npoint,chirplen,&maxi,&snr_max,&lin0,&lin90,&var);

      /* identify when an impulse would have caused observed filter output */
      impulseoff=(maxi+chirppoints)%npoint;
      timeoff=datastart+impulseoff/srate;
      timestart=datastart+maxi/srate;

      /* if SNR greater than 5, then print details, else just short message */
      if (snr_max<5.0)
         printf("max snr: %.2f offset: %d data start: %.2f sec. variance: %.5f\n",
		snr_max,maxi,datastart,var);
      else {
         /* See if the nominal chirp can pass a frequency-space veto test */
         varsplit=splitup_freq2(lin0*n0/sqrt(2.0),lin90*n90/sqrt(2.0),ch0tilde,ch90tilde,2.0/(n0*n0),
                                twice_inv_noise,npoint,maxi,8,indices,stats,output0,htilde);
         prob=gammq(4.0,4.0*varsplit);
         printf("\nMax SNR: %.2f (offset %d) variance %f\n",snr_max,maxi,var);
         printf("   If impulsive event, offset %d or time %.2f\n",impulseoff,timeoff);
         printf("   If inspiral, template start offset %d (time %.2f) coalescence time %.2f\n",
                         maxi,timestart,timestart+tc);
         printf("   Normalization: S/N=1 at %.2f kpc\n",1000.0*distance);
         printf("   Linear combination of max SNR: %.4f x phase_0 + %.4f x phase_pi/2\n",lin0,lin90);
         if (prob<0.01)
            printf("   Less than 1%% probability that this is a chirp (p=%f).\n",prob);
         else
            printf("   POSSIBLE CHIRP!  with > 1%% probability (p=%f).\n",prob);

         /* See if the time-domain statistics are unusual or appears Gaussian */
         if (is_gaussian(datas,npoint,-2048,2047,1))
            printf("   Distribution does not appear to have outliers...\n\n");
         else
            printf("   Distribution has outliers! Reject\n\n");
      }
   
      /* shift ends of buffer to the start */
      needed=npoint-chirplen+1;
      for (i=0;i<chirplen-1;i++) datas[i]=datas[i+needed];
   }
   return 0;
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
      fgetinput.npoint=fgetoutput.srate*min_into_lock*60.0-needed;
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

   num_sent++;
   if (reset==1)
      return 1;
   else
      return 2;
}
