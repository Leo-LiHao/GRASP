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
  static FILE *fpifo,*fplock;
  static int first=1,remain=0,num_sent=0;
  float tstart;
  int i,temp,code=2,diff=NPOINT-needed;

  if (first) { /* on first call only */
    FILE *fpss;
    first = 0;
    diff = 0;
    /* open the IFO output file, lock file, and swept-sine file */
    fpifo = grasp_open("GRASP_DATAPATH","channel.0","r");
    fplock = grasp_open("GRASP_DATAPATH","channel.10","r");
    fpss = grasp_open("GRASP_DATAPATH","swept-sine.ascii","r");
    /* get the response function and put in scaling factor */
    normalize_gw(fpss,NPOINT,srate,response);
    for (i=0;i<NPOINT;i++) response[i] *= HSCALE/ARMLENGTH_1994;
    fclose(fpss);
  }

  if (num_sent==DATA_SEGMENTS) return 0;

  /* if new locked section, skip forward */
  while (remain<needed) {
    fprintf(stderr,"\nEntering new locked set of data\n");
    temp = get_data(fpifo,fplock,&tstart,MIN_INTO_LOCK*60*srate,datas,&remain,&srate,1);
    if (temp==0) return 0;

    /* number of points needed will be full length */
    needed = NPOINT;
    diff = 0;
    code = 1;
  }

  /* get the needed data and compute the start time of the buffer */
  temp = get_data(fpifo,fplock,&tstart,needed,datas+diff,&remain,&srate,0);
  if (temp==0) return 0;
  datastart = tstart - diff/srate;
  
  num_sent++;
  return code;
}
