/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"

#define NPOINT 131072      /* number of data points */
#define HSCALE 1.0e21      /* convenient scaling factor */
#define FLO 120.0          /* low frequency cutoff for filtering */
#define MIN_INTO_LOCK 3.0  /* time (minutes) to skip into each locked section */
#define THRESHOLD 6.0      /* detection threshold SNR */
#define ATTEN 30.0         /* attenuation cutoff for ringdown waveforms */
#define SAFETY 1000        /* padding safety to avoid wraparound errors */
#define DATA_SEGMENTS 3000 /* maximum number of data segments to filter */

double datastart;
float *response,srate=9868.4208984375;
short *datas;
int needed=NPOINT;

int main()
{
  void realft(float *, unsigned long, int);
  int fill_buffer();
  double norm;
  float *data,*htilde,*output;
  float *mean_pow_spec,*twice_inv_noise;
  float *ring,*ringtilde,*template;
  float decaytime,decay=0.0,scale,snr,mean,var,tmpl_norm,dist;
  float mass=50.0,spin=0.98,eps=0.03,psi0=0.0,invMpc=10.0,ringstart=500.0;
  int i,code,len,safe=SAFETY,diff,off,n=NPOINT;

  /* allocate memory for arrays */
  response=(float *)malloc(sizeof(float)*(NPOINT+2));
  datas=(short *)malloc(sizeof(short)*NPOINT);
  data=(float *)malloc(sizeof(float)*NPOINT);
  htilde=(float *)malloc(sizeof(float)*NPOINT);
  output=(float *)malloc(sizeof(float)*NPOINT);
  ringtilde=(float *)malloc(sizeof(float)*NPOINT);
  template=(float *)malloc(sizeof(float)*NPOINT);
  mean_pow_spec=(float *)malloc(sizeof(float)*(NPOINT/2+1));
  twice_inv_noise=(float *)malloc(sizeof(float)*(NPOINT/2+1));
  
  /* manufacture quasinormal ring data; obtain length of signal */
  ring = NULL;
  len = qn_qring(psi0,eps,mass,spin,1.0/srate,ATTEN,n,&ring);

  /* normalize quasinormal ring to one megaparsec */
  scale = HSCALE*M_SOLAR/MPC;
  for (i=0;i<len;i++) ringtilde[i] = ring[i] *= scale;
  for (i=len;i<n;i++) ringtilde[i] = ring[i] = 0;

  /* FFT the quasinormal ring waveform */
  realft(ringtilde-1,n,1);
  if (n<len+2*safe) abort();

  while (1) {

    /* fill buffer with number of points needed */
    code = fill_buffer();

    /* if no points left, we are done! */
    if (code==0) break;

    /* if just entering a new locked stretch, reset averaging over power spectrum */
    if (code==1) {
      norm = 0;
      clear(mean_pow_spec,n/2+1,1);

      /* decay time in seconds: set to 15 x length of NPOINT sample */
      decaytime = 15.0*n/srate;
      decay = exp(-1.0*n/(srate*decaytime));
    }

    /* copy data into floats */
    for (i=0;i<NPOINT;i++) data[i] = datas[i];

    /* inject a time-domain signal before FFT (note output is used as temp storage only) */
    qn_inject(data,ring,response,output,invMpc,(int)(srate*(ringstart-datastart)),n,len);

    /* compute the FFT of data */
    realft(data-1,n,1);

    /* normalized dL/L tilde */
    product(htilde,data,response,n/2);

    /* update auto-regressive mean power spectrum */
    avg_inv_spec(FLO,srate,n,decay,&norm,htilde,mean_pow_spec,twice_inv_noise);

    /* normalize the ring to produce a template */
    qn_normalize(template,ringtilde,twice_inv_noise,n,&tmpl_norm);

    /* calculate the filter output and find its maximum */
    find_ring(htilde,template,twice_inv_noise,output,n,len,safe,&off,&snr,&mean,&var);

    /* perform diagnostics on filter output */
    if (snr<THRESHOLD) { /* threshold not exceeded: print a short message */
      printf("max snr: %.2f (offset %6d) ",snr,off);
      printf("data start: %.2f variance: %.5f\n",datastart,var);
    } else { /* threshold exceeded */
      /* estimate distance to signal (template distance [Mpc] = 1 / tmpl_norm) */
      dist = 2/(tmpl_norm*snr);
      printf("\nMax SNR: %.2f (offset %d) variance %f\n",snr,off,var);
      printf("   If ringdown, estimated distance: %f Mpc, ",dist);
      printf("start time: %f\n",datastart+off/srate);
      /* See if time domain statistics are non-Gaussian */
      if (is_gaussian(datas,n,-2048,2047,1))
	printf("   POSSIBLE RINGDOWN: Distribution does not appear to have outliers\n\n");
      else
	printf("   Distribution has outliers!  Reject\n\n");
    }

    /* shift ends of buffer to the start */
    diff = len + 2*safe; /* safety is applied at beginning and end of buffer */
    needed = NPOINT - diff;
    for (i=0;i<diff;i++) datas[i] = datas[i+needed];
  }

  return 0;
}

/* this routine gets the data, overlapping the data buffer as needed */
int fill_buffer()
{
  static struct fgetinput fgetinput;
  static struct fgetoutput fgetoutput;
  static double lastcalibtime=0;
  static int first=1,num_sent=0;
  int i,temp,code=2,diff=NPOINT-needed;

  if (first) { /* on first call only */
    first = 0;
    /* stores ADC data as short integers */
    fgetinput.nchan = 1;
    fgetinput.files = framefiles;
    fgetinput.calibrate = 1;
    fgetinput.chnames = (char **)malloc(fgetinput.nchan*sizeof(char *));
    fgetinput.locations = (short **)malloc(fgetinput.nchan*sizeof(short *));
    fgetoutput.npoint = (int *)malloc(fgetinput.nchan*sizeof(int));
    fgetoutput.ratios = (int *)malloc(fgetinput.nchan*sizeof(int));
    fgetoutput.lastlock = fgetoutput.tstart = 0;
    /* set up channel names */
    fgetinput.chnames[0] = "IFO_DMRO";
    if (NULL!=getenv("GRASP_REALTIME")) fgetinput.inlock = 0; /* don't care if locked */
    else fgetinput.inlock = 1; /* only locked */
  }

  if (num_sent==DATA_SEGMENTS) return 0;

  do {

    /* get number of points required */
    fgetinput.npoint = needed;
    fgetinput.locations[0] = datas + diff;
    fgetinput.seek = 0;

    temp = fget_ch(&fgetoutput,&fgetinput);
    srate = fgetoutput.srate;
    datastart = fgetoutput.dt - diff/srate;

    /* if nothing left, return */
    if (temp==0) return 0;

    /* if new locked section, skip forward */
    if (temp==1) {
      fprintf(stderr,"\nEntering new locked set of data\n");
      fgetinput.seek = 1;
      fgetinput.npoint = MIN_INTO_LOCK*60*srate - needed;
      if (0==fget_ch(&fgetoutput,&fgetinput)) return 0;

      /* number of points needed will be full length */
      needed = NPOINT;
      diff = 0;
      code = 1;
    }

  } while (temp==1);

  /* if we have recalibrated, get the response function and put in scaling factor */
  if (lastcalibtime<fgetoutput.tcalibrate) {
    lastcalibtime = fgetoutput.tcalibrate;
    GRnormalize(fgetoutput.fri,fgetoutput.frinum,NPOINT,srate,response);
    for (i=0;i<NPOINT;i++) response[i] *= HSCALE/ARMLENGTH_1994;
  }
  
  num_sent++;
  return code;
}
