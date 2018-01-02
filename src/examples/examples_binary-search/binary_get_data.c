#include <string.h>
#include "grasp.h"
#include "binary_params.h"

#define SIM_VARIANCE 16384.0
#define SIM_SITE 8

char *rcsid="$id";

static int count_chunks,count_locked,count_segments;
static double count_tinlock,lastlock,discarded,lastdiscard=0;
static FILE *fp_statistics;

struct Chunk {
    double time;
    short *data;
    float *spec;
    int cont;
    int is_gaussian;
    int counter;
    int used_in_spectrum;
};

void make_retlifs(float m1, float m2, float *ch1, float *ch2, float fstart,
                  int n, float srate, int *filled, float *t_coal,
                  int err_cd_sprs, int order);

int get_calibrated_data();
#if (RANDOMIZE)
float ran2(long *);
long randomize=-12345;
#endif

void nullout(float freqmin,float freqmax,int npoint,float srate,float* n_inv_noise) {
  int imin,imax,i;
  
  imin=freqmin*npoint/srate;
  imax=freqmax*npoint/srate;
  /* set to zero in those bins */
  if (imin<npoint/2)
    for (i=imin;((i<imax) && (i<npoint/2));i++)
      n_inv_noise[i]=0.0;
  return;
}

/* Routine to compute window function for some specified window type.
 * 
 * Input is the window type:
 *   type = 0:  rectangular (no) window,
 *   type = 1:  Hann window,
 *   type = 2:  Welch window,
 *   type = 3:  Bartlett window.
 *   
 * Output are *wss (sum of window function values squared) and the window
 * function values window[0..n-1]. */
void compute_window(float *wss, float window[], int n, int type)
{
  int i;
  float fac=2.0/n;
  
  *wss = 0;
  for (i=0;i<n;i++) {
    float win;
    switch (type) {
    case 0:    /* rectangular (no) window */
      win = 1;
      break;
    case 1:    /* Hann window */
      win = 0.5*(1 - cos(i*fac*M_PI));
      break;
    case 2:    /* Welch window */
      win = i*fac - 1;
      win = 1 - win*win;
      break;
    case 3:    /* Bartlett window */
      win= 1 - fabs(i*fac - 1);
      break;
    default:   /* unrecognized window type */
      GR_start_error("avg_spec()",rcsid,__FILE__,__LINE__);
      GR_report_error("don't recognize windowtype=%d\n",type);
      GR_end_error();
      abort();
      break;
    }
    *wss += win*win;
    window[i] = win;
  }
  return;
}

/* Routine to calculate the response function using the old-style
 * data acquisition routines */
void recalibrate_old(float *response, int npoint)
{
  int i;
  FILE *fpss;
  fpss = grasp_open("GRASP_DATAPATH","swept-sine.ascii","r");
  normalize_gw(fpss,npoint,SRATE,response);
  fclose(fpss);
  for (i=0;i<npoint;i++) response[i] *= HSCALE/ARMLENGTH_1994;
  return;
}


/* Are we building a frame-compatible version? */
#if (DATA_FORMAT == 1)

/* Global frame variables */
struct fgetoutput fgetoutput={0};
struct fgetinput fgetinput={0};
int frame_init=0;
int time_last_calibrated=-1;

/* Routine to initialize global frame variables and to allocate memory */
void initialize_frame()
{
  frame_init = 1;
  fgetinput.nchan = 1;
  fgetinput.files = framefiles;
  fgetinput.calibrate = 1;
  fgetinput.chnames = (char **)malloc(fgetinput.nchan*sizeof(char *));
  fgetinput.locations = (short **)malloc(fgetinput.nchan*sizeof(short *));
  fgetinput.chnames[0] = "IFO_DMRO";
  fgetinput.inlock = 1;
  fgetoutput.npoint = (int *)malloc(fgetinput.nchan*sizeof(int));
  fgetoutput.ratios = (int *)malloc(fgetinput.nchan*sizeof(int));
  return;
}

/* Routine to calculate the response function using the frame
 * data acquisition routines.  This routine should only be called
 * after some data has been read, so that the frame output variables
 * have been set. */
void recalibrate_frame(float *response, int npoint)
{
  float srate=SRATE;
  int i;
  
  if (!frame_init) {
    GR_start_error("recalibrate_frame()",rcsid,__FILE__,__LINE__);
    GR_report_error("must get some data before calling this routine!\n");
    GR_end_error();
    abort();
  }
  
  /* print message */
  GR_start_error("recalibrate_frame()",rcsid,__FILE__,__LINE__);
  GR_report_error("Note: data calibration carried out at time %d.\n",
		  fgetoutput.tcalibrate);
  GR_report_error("  previous frame calibration was from time %d.\n",
		  time_last_calibrated);
  GR_end_error();

  /* do the calibration */
  GRnormalize(fgetoutput.fri,fgetoutput.frinum,npoint,srate,response);
  for (i=0;i<npoint;i++) response[i] *= HSCALE/ARMLENGTH_1994;

  /* record time that we did calibration */
  time_last_calibrated=fgetoutput.tcalibrate;
  
  return;
}
#endif  /* end of frame-compatible version */


/* Routine to calculate the response function for simulated data */
void recalibrate_sim(float *response, int npoint)
{
  double *power,*tmp;
  float parameters[9],srate=SRATE;
  float fac=HSCALE*sqrt(0.5*srate/SIM_VARIANCE);
  char  site_name[256],noise_file[256],whiten_file[256];
  int m=npoint/2;
  tmp = power = (double *)malloc(m*sizeof(double));
  if (power == NULL) {
    GR_start_error("recalibrate_sim()",rcsid,__FILE__,__LINE__);
    GR_report_error("failed to allocate %d doubles\n",m);
    GR_end_error();
    abort();
  }
  detector_site("detectors.dat",SIM_SITE,parameters,
		site_name,noise_file,whiten_file);
  noise_power(noise_file,m,srate/npoint,power);
  while (m-- > 0) {
    *response++ = fac*sqrt(*tmp++);
    *response++ = 0;
  }
  free(power);
}


/* Routine to calculate the response function.  The response function
 * is only updated if needed --- otherwise, it is unchanged.
 * 
 * Input is the number of points in the response function, npoint.
 * 
 * Output is the response function array response[0..npoint-1]. */
void recalibrate(float *response, int npoint)
{
# if (DATA_FORMAT == 0)
    static int first = 1;
    /* old format data */
    if (first) {
      first = 0;
      recalibrate_old(response,npoint);
    }
# elif (DATA_FORMAT == 1)
    /* frame format data */
    if (time_last_calibrated != fgetoutput.tcalibrate)
      recalibrate_frame(response,npoint);
# elif (DATA_FORMAT == 2)
    static int first = 1;
    /* simulated data */
    if (first) {
      first = 0;
      recalibrate_sim(response,npoint);
    }
# else
    /* not a recognized format! */
    GR_start_error("recalibrate()",rcsid,__FILE__,__LINE__);
    GR_report_error("unrecognized DATA_FORMAT %d\n",DATA_FORMAT);
    GR_end_error();
    abort();
# endif
  return;
}

/* open the IFO output file and read the mainheader to get runstart time */
double get_runstart_old()
{
  struct ld_mainheader mh;
  struct ld_binheader bh;
  float tstart,srate;
  int zero=0,npt=NPOINT;
  FILE *fpifo;
  short *here;
  fpifo = grasp_open("GRASP_DATAPATH","channel.0","r");
  read_block(fpifo,&here,&zero,&tstart,&srate,0,&npt,1,&bh,&mh);
  fclose(fpifo);
  return (double)mh.epoch_time_sec + (double)mh.epoch_time_msec*0.001;
}

/* Routine to get a chunk of data using the old data acquisition routines.
 * 
 * Input is the number of points in the chunk, npoint.
 * 
 * Output is the data in the array here[0..npoint-1]
 * and the time, *time, of the start of the data.
 * 
 * Returned is the code:
 *   code = 0: no more data;
 *   code = 1: beginning of a locked section
 *             (no data acquired, but we have skipped MIN_INTO_LOCK
 *              minutes into the locked segment);
 *   code = 2: continuing a locked section (data acquired). */
int get_chunk_data_old(short *here, int npoint, double *time)
{
  static FILE *fpifo,*fplock;
  static int remain=0,first=1;
  static float srate=SRATE;
  static double runstart;
  float tstart;
  int code;
  
  if (first) { /* initialization */
    first = 0;
    runstart = get_runstart_old();
    /* open the IFO output file, lock file, and swept-sine file */
    fpifo = grasp_open("GRASP_DATAPATH","channel.0","r");
    fplock = grasp_open("GRASP_DATAPATH","channel.10","r");
  }

  if (remain < npoint) { /* if new locked section */
    int nskip=(int)(60*MIN_INTO_LOCK*SRATE);
    /* skip forward MIN_INTO_LOCK minutes */
    code = get_data(fpifo,fplock,&tstart,nskip,here,&remain,&srate,1);
    if (code == 0) return 0;
    else return 1;
  }

  /* get the needed data */
  code = get_data(fpifo,fplock,&tstart,npoint,here,&remain,&srate,0);
  
  /* compute the start time */
  *time = runstart + tstart;

  return code;
}


/* are we building a frame compatible version? */
#if (DATA_FORMAT == 1)

/* Routine to get a chunk of data using the frame data acquisition routines.
 * 
 * Input is the number of points in the chunk.
 * 
 * Output is the data in the array here[0..npoint-1]
 * and the time, *time, of the start of the data.
 * 
 * Returned is the code:
 *   code = 0: no more data;
 *   code = 1: beginning of a locked section
 *             (no data acquired, but we have skipped MIN_INTO_LOCK
 *              minutes into the locked segment---WARNING... the data
 *              array here will be filled with junk, which should be
 *              ignored!)
 *   code = 2: continuing a locked section (data acquired). */
int get_chunk_data_frame(short *here, int npoint, double *time, int *counter)
{
  static float srate=SRATE;
  int code, retcode;
  static int private_count=0;

  private_count++;
  fprintf(fp_statistics,"get_chunk_data_frame: getting chunk %d\n",private_count);
  
  if (!frame_init) initialize_frame();
  
                                   /* we want... */
  fgetinput.npoint = npoint;       /* ...npoints of data... */
  fgetinput.locations[0] = here;   /* ...in array here[]... */
  fgetinput.seek = 0;              /* ...no seek... */
  fgetinput.inlock = 1;            /* ...only when in lock */

  /* get npoints of data */
  retcode = code = fget_ch(&fgetoutput,&fgetinput);
  fprintf(fp_statistics,"Current frame file: %s\n",fgetoutput.filename);

  while (code == 1) {
    /* seek to a total of MIN_INTO_LOCK minutes into locked section ... */
    fgetinput.seek = 1;
    /* ... we have already gone npoint, so we need to subtract this off ... */
    fgetinput.npoint = MIN_INTO_LOCK*60*srate - npoint;
    /* ... subtracting npoint is an inconvenience! */
    if (fgetinput.npoint < 1) {
      GR_start_error("get_chunk_data_frame()",rcsid,__FILE__,__LINE__);
      GR_report_error("set MIN_INTO_LOCK to a larger value\n");
      GR_end_error();
      abort();
    }
    code = fget_ch(&fgetoutput,&fgetinput);
  fprintf(fp_statistics,"Current frame file: %s\n",fgetoutput.filename);

    if (code == 0) return 0; /* no more data */
    if (code == 1) { /* oops... skiped into yet another locked section ... */
      /* ... so we need to skip forward the npoint offset we subtracted above */
      fgetinput.npoint = npoint;
      code = fget_ch(&fgetoutput,&fgetinput);
  fprintf(fp_statistics,"Current frame file: %s\n",fgetoutput.filename);

      if (code == 0) return 0;
    }
  }

  lastlock=fgetoutput.lastlock;
  discarded=(double)fgetoutput.discarded/srate;
  *time = fgetoutput.tstart;
  *counter=private_count;
  return retcode;
}
#endif

/* Routine to generate a chunk of simulated white noise.
 * 
 * Input is the number of points in the chunk, npoint.
 * 
 * Output is the data in the array here[0..npoint-1]
 * and the time, *time, of the start of the data.
 * 
 * Returned is the code:
 *   code = 0: no more data;
 *   code = 1: beginning of a locked section
 *             (no data acquired, but we have skipped MIN_INTO_LOCK
 *              minutes into the locked segment);
 *   code = 2: continuing a locked section (data acquired). */
int get_chunk_data_sim(short *here, int npoint, double *time)
{
  float ran2(long *);
  static double local_time=0,srate=SRATE;
  static long seed=-100;
  int m=npoint/2;
  
  if (seed < 0) {
    ran2(&seed);
    return 1;
  }
  local_time += (double)npoint/srate;
  *time = local_time;
  while (m-- > 0) {
    float x,y,rr,fac;
    do {
      x = 2*ran2(&seed) - 1;
      y = 2*ran2(&seed) - 1;
      rr = x*x + y*y;
    } while (rr > 1 || rr == 0);
    fac = sqrt(-2*SIM_VARIANCE*log(rr)/rr);
    *here++ = floor(fac*x);
    *here++ = floor(fac*y);
  }
  return 2;
}


/* Routine to compute the spectrum of a chunk of data.
 * 
 * Input is the number of points of data, npoint,
 * and the array of data, data[0..npoint-1].
 * 
 * Output is spectrum, spec[0..npoint/2]. */
void get_spectrum(short *data, float *spec, int npoint)
{
  void realft(float [], unsigned long, int);
  static float *work=NULL,*win=NULL,fac;
  static int nkeep=-1;
  int i,n=npoint,m=npoint/2;
  
  if (npoint != nkeep) {
    float wss;
    nkeep = npoint;
    work = realloc(work,npoint*sizeof(float));
    if (work == NULL) {
      GR_start_error("get_spectrum()",rcsid,__FILE__,__LINE__);
      GR_report_error("failed to allocate %d floats\n",npoint);
      GR_end_error();
      abort();
    }
    win = realloc(win,npoint*sizeof(float));
    if (win == NULL) {
      GR_start_error("get_spectrum()",rcsid,__FILE__,__LINE__);
      GR_report_error("failed to allocate %d floats\n",npoint);
      GR_end_error();
      abort();
    }
    compute_window(&wss,win,npoint,WINDOW);
    fac = 2/(wss*SRATE);
  }
  /* copy data to workspace and apply window */
  for (i=0;i<n;i++) work[i] = win[i]*data[i];
  
  /* FFT workspace */
  realft(work-1,npoint,1);
  
  /* DC and Nyquist terms */
  spec[0] = fac*work[0];
  spec[m] = fac*work[1];
  
  /* other terms */
  for (i=1;i<m;i++) {
    int ir=i+i,ii=ir+1;
    float re=work[ir],im=work[ii];
    spec[i] = fac*(re*re + im*im);
  }
  return;
}

/* === Patrick Brady === Modified: 4 December 1998 === 
 *
 * Routine to fill a chunk with data.
 * 
 * Input is the number of points, npoint, in the chunk of data.
 * 
 * Output is the struct chunk with
 *   chunk.data[0..npoint-1]: the data
 *   chunk.spec[0..npoint/2]: the power spectrum of the data
 *   chunk.time: the time of the start of the data in the chunk
 *   chunk.cont: a flag indicating if the data is continuous from the
 *     last chunk filled (1) or not (0).
 *   chunk.is_gaussian: a flag to indicate if the data has outliers (0) or not (1)
 * 
 * Returned is the code:
 *   code = 0: no more data;
 *   code = 1: beginning of a locked section
 *             (no data acquired, but we have skipped MIN_INTO_LOCK
 *              minutes into the locked segment);
 *   code = 2: continuing a locked section (data acquired). */
int fill_chunk(struct Chunk *chunk, int npoint)
{
  int code;
  
#if (DATA_FORMAT == 0)
  /* old format data */
  code = get_chunk_data_old(chunk->data,npoint,&chunk->time);
#elif (DATA_FORMAT == 1)
  /* frame format data */
  code = get_chunk_data_frame(chunk->data,npoint,&chunk->time,&chunk->counter);
#elif (DATA_FORMAT == 2)
  /* simulated data */
  code = get_chunk_data_sim(chunk->data,npoint,&chunk->time);
#else
  /* not a recognized request! */
  GR_start_error("fill_chunk()",rcsid,__FILE__,__LINE__);
  GR_report_error("unrecognized DATA_FORMAT %d\n",DATA_FORMAT);
  GR_end_error();
  abort();
#endif

  switch (code) {
  case 0: /* no more data */
      return code;
  case 1: /* starting a new locked section */
      chunk->cont = 0;
      break;
  case 2: /* continuing a locked section */
      chunk->cont = 1;
      break;
  default:
      GR_start_error("fill_chunk()",rcsid,__FILE__,__LINE__);
      GR_report_error("unrecognized code %d\n",code);
      GR_end_error();
      abort();
  }

  chunk->is_gaussian=is_gaussian(chunk->data,npoint,-2048,2047,0);
  chunk->used_in_spectrum=0; 
  get_spectrum(chunk->data,chunk->spec,npoint);
  
  return code;
}

/* Bruce Allen, Modified Jan 18, 1998
   The new algorithm works as follows:
   (1) go through all the kk current chunks
   (2) Count the number which are outlier free
   (3) replace the oldest of the kk saved spectra with these ones
   (4) compute the average spectrum
   (5) return

   Comment: the only place this routine might "misbehave" is at the start of
   a newly-calibrated section of data (if the calibration curve has significantly
   changed), or if the first kk segments of data all have
   outliers.
*/
void average_spectrum(struct Chunk *chunk, int kk, int npoint, float *spec, int new_lock)
{
    static int num_stored=0,oldest=0;
    static float *stored_spectras=NULL;
    static int *stored_chunks=NULL;
    
    float factor;
    int i,j,m=npoint/2;
    
        /* allocate memory if first time */
    if (!stored_spectras) 
        {
        stored_spectras=(float *)malloc(kk*m*sizeof(float));
        stored_chunks=(int *)malloc(kk*sizeof(int));
        }
    
        /* loop over outlier-free chunks not already used in spectrum */
    for (i=0;i<kk;i++)
        if (chunk[i].is_gaussian && !chunk[i].used_in_spectrum) 
            {
                float *oldest_spec;
                    /* mark this spectrum as used */
                chunk[i].used_in_spectrum=1;
                    /* keep track of chunk being used to produce spectrum */
                stored_chunks[oldest]=chunk[i].counter;
                    /* copy outlier-free spectra to oldest saved spectra */
                oldest_spec=stored_spectras+oldest*m;
                memcpy((void *)(oldest_spec),(const void *)(chunk[i].spec),(size_t)m*sizeof(float));
                    /* make circular buffer point back if needed */
                oldest=(oldest+1)%kk;
                    /* increment number of stored spectra */
                if (num_stored<kk) num_stored++;
                    /* print some useful diagnostic information */
                fprintf(fp_statistics,"Incorporating chunk %d into power spectra",chunk[i].counter);
                for (j=0;j<num_stored-1;j++)
                    fprintf(fp_statistics," %d",stored_chunks[j]);
                fprintf(fp_statistics," %d\n",stored_chunks[num_stored-1]);
            }
    
        /* check that we have SOMETHING stored! */
    if (!num_stored)
        {        
            fprintf(stderr,"average_spectrum: no outlier free data to compute spectrum with!\n");
            fflush(stderr);
            abort();
        }
    
        /* Now compute the sum of the stored spectra */
    clear(spec,m,1);
    for (j=0;j<num_stored;j++) 
        {
            float *current_spec;
            current_spec=stored_spectras+j*m;
            for (i=0;i<m;i++)
                spec[i]+=current_spec[i];
        }
    
        /* Normalize to get average spectrum  */
    factor=1.0/num_stored;
    for (i=0;i<m;i++)
        spec[i]*=factor;

        /* print out some information about which chunks are used */
   
    return;
}

/* === Patrick Brady === Modified: 4 December 1998 ===
  
   Routine to average the spectra of several chunks of data.
   
   chunk: Input. The array of chunks

   kk: Input. The number of chunks

   npoint: Input. The number of points of data in each chunk

   spec: Ouput. Pointer to the array of (npoint/2 +1) floats

   new_lock: Input. Flag to indicate if all the chunks are continuous (new_lock=0)
             or if we have new locked data in the buffer (new_lock=1).

   Authors: Jolien Creighton
  
 */
void average_spectrum_brady(struct Chunk *chunk, int kk, int npoint, float *spec, int new_lock)
{
    static float *stored_spec;
    static int first=1;
    float fac=1/(float)kk,norm=0;
    int i,j,m=npoint/2;
    
    if (first){ /* Allocate memory the first time in */
        stored_spec=(float *)malloc((npoint/2 + 1)*sizeof(float));
        clear(stored_spec,m,1);
        first=0;
    }
        /* Note: ignore the Nyquist component for some reason.... */
        /* If the data is continuous ... */
    if (new_lock==0){
        clear(spec,m,1);          
        for (j=0;j<kk;j++){
                /* If data is outlier free ... */
            if (chunk[j].is_gaussian){
                    /* ..... add in spectrum from this chunk */
                norm+=1;
                for (i=0;i<m;i++)
                    spec[i] += chunk[j].spec[i];  
            }
        }
            /* If all spectra have outliers ... */
        if (!norm){
            GR_start_error("average_spectrum()",rcsid,__FILE__,__LINE__);
            GR_report_error("%i contiguous chunks with outliers\n",kk);
            GR_end_error();
                /* return the spectrum computed previously */
            for (i=0;i<m;i++) spec[i] = stored_spec[i]; 
        }
            /* otherwise ... */
        else{
            fac= (1/norm);
            for (i=0;i<m;i++){
                spec[i] *= fac; 	 /* normalize the spectrum */
                stored_spec[i]= spec[i]; /* keep this spectrum for when we encounter new locked data */
            }
        }
    }
        /* otherwise the data is not continuous */
    else 
        for (i=0;i<m;i++)
            spec[i] = stored_spec[i]; /* return the spectrum computed previously */
    return;
}

/* === Patrick Brady === 4 December 1998 ===

   This function initialises the array of chunks,  and returns an
   integer code:
      code = 0: no more data;
      code = 2: continuing a locked section.

   If it returns code = 1 something is wrong.
   
   The arguments are:
   
     chunk: Output. The array of chunks.

     npoint: Input. The number of points in chunk.data[0..npoint-1]

     chunks_filled: Input. Number of chunks that are already filled.

     start_filled: Input. First chunk with good data

     kk: Input: total number of chunks

*/

int initialise_chunks(struct Chunk *chunk, int npoint, int chunks_filled, int start_filled, int kk){
    int i=0,code=2;
    float srate=SRATE;

    while (i < chunks_filled){
	/* copy the data */
	memcpy((void *)chunk[i].data,(void *)chunk[start_filled+i].data,(size_t)npoint*sizeof(short));	
	/* copy the spectrum */
	memcpy((void *)chunk[i].spec,(void *)chunk[start_filled+i].spec,(size_t)(npoint/2+1)*sizeof(float));
	/* and the time,  and continuity flag */
	chunk[i].time=chunk[start_filled+i].time;
	chunk[i].cont=chunk[start_filled+i].cont;
    chunk[i].counter=chunk[start_filled+i].counter;
    chunk[i].is_gaussian=chunk[start_filled+i].is_gaussian;
    chunk[i].used_in_spectrum=chunk[start_filled+i].used_in_spectrum;
	/* increment the counter */
	i++;
    }

    while (i < kk){
	code = fill_chunk(chunk + i,npoint);
	switch (code) {
	case 0:  /* no more data */
	    return 0;
	case 1:  /* entering new locked set */
	    fprintf(fp_statistics,"Message from initialise_chunks():\n");
	    fprintf(fp_statistics,"  Locked section was too short to get %i chunks\n",kk);
	    fflush(fp_statistics);
	    i= -1;
	    break;
	case 2:  /* continuing a locked set */
	    break;
	default: /* unrecognized code */
	    GR_start_error("initialise_chunks()",rcsid,__FILE__,__LINE__);
	    GR_report_error("unrecognized code %d\n",code);
	    GR_end_error();
	    abort();
	}
	/* increment the counter */
	i++;
    }

    /* If we have a locked section > MIN_INTO_LOCK + kk*npoint/srate */
    if(code==2){ 
	fprintf(fp_statistics,"Message from initialise_chunks():\n");
	fprintf(fp_statistics,"  Last locked section had %f secs of data, ",MIN_INTO_LOCK*60.0 +
		(float)(count_chunks)*npoint/srate);
	fprintf(fp_statistics," and %f secs discarded\n",discarded-lastdiscard);
	lastdiscard=discarded;
	count_locked++;
	count_chunks=kk;
	count_tinlock+=(MIN_INTO_LOCK*60.0 + (double)(kk*npoint)/srate);
	fprintf(fp_statistics,"  Starting locked section %i at time %f with segment %d\n",count_locked,lastlock,count_segments);
	fflush(fp_statistics); 
    } 

    return code;
}

/* === Patrick Brady === Modified: 4 December 1998 ===

   This function gets more data for the binary_search code.   It returns an
   integer code:
      code = 0: no more data;
      code = 1: beginning of (well... MIN_INTO_LOCK minutes into)
                a locked section of data
      code = 2: continuing a locked section.
   
   The arguments are:
   
     need: Input.   The number of points to skip ahead between calls.

     data: Ouput.   The IFO_DMRO is returned in this array

     spec: Output.  The average spectrum is returned in spec[0..npoint/2]

     *time: Output. Time at the start of the data

   Authors:  Jolien Creighton.
   
*/

int get_more_data(int need, short data[], float spec[], float response[],
                  double *time)
{
  static short *buffer;              /* pointer to the beginning of the buffered data      */
  static short *end_buffer;          /* pointer to the end of the buffered data            */
  static short *here;                /* pointer to the position of the data to return next */
  static struct Chunk *chunk;        /* buffered data also stored in chunks                */
  static int lastfill_chunkNumber=-1;/* last chunk that was filled with data               */
  static int first=1,npoint=NPOINT,chunk_needed=0,here_position=0,new_lock=0;
  const int kk=8,k=4;                /* kk is the (even) number of chunks; k is half of kk */
  int here_chunkNumber=0,out_code=2; /* chunk that here points to,  and the return code    */
  
  if (first) {/* on first call */
    int i;
    first = 0;
    /* Open the file to record the statistics of the run */
    fp_statistics = grasp_open("GRASP_MFPATH","run.stats","w");
    count_chunks=count_locked=count_segments=0;
    count_tinlock=0.0;
    /* allocate memory to the buffer: extra npoints is to insure continuity in the data */
    here = buffer = (short *)malloc((kk + 1)*npoint*sizeof(short));
    /* end_buffer points to the beginning of the duplicated data at the end of the buffer */
    end_buffer = buffer + kk*npoint;
    /* allocate memory for the chunks of data */
    chunk = (struct Chunk *)malloc(kk*sizeof(struct Chunk));
    for (i=0;i<kk;i++) {  /* partition the buffered data into chunks */
      chunk[i].data = buffer + i*npoint;
      chunk[i].spec = (float *)malloc((npoint/2 + 1)*sizeof(float));
    }
    /* fill the chunks with data */
    if(!(initialise_chunks(chunk,npoint,0,0,kk))) return 0;
  } 
  else {/* on other calls */
      here += need;            /* advance here by the needed amount of data       */
      chunk_needed += need;    /* increment flag for getting another chunk        */
      here_position += need;   /* increment counter to position of here in buffer */
      here_chunkNumber = (here_position/npoint)%kk; /* chunk that here points to  */
  }
  
  /* If here is in duplicated data,  move to equivalent place at start of buffer */
  if ( (here_position/npoint) == kk ){
      here_position = here_position%npoint;
      here = buffer + here_position;
  }

  /* Does next chunk contain continuous data? */
  if(here_position%npoint > 0){
      int next_chunk=(here_chunkNumber + 1)%kk;
      /* If not ..... */
      if (chunk[next_chunk].cont == 0 && new_lock == 1){
	  int code;
	  code = initialise_chunks(chunk,npoint,0,0,kk);
	  switch (code) {
	  case 0:  /* no more data */
	      return 0;
	  case 2:  /* continuing a locked set */
	      break;
	  default: /* unrecognized code */
	      GR_start_error("get_more_data()",rcsid,__FILE__,__LINE__);
	      GR_report_error("Should never get code %d here\n",code);
	      GR_end_error();
	      abort();
	  }
	  /* PATRICK:  check this with Jolien and Bruce.  Does the calling routine need 
	     to know about this?   Or should it be told code=2?
	  */
	  out_code=1;                /* let calling routine know this is new locked data */
	  here=buffer;               /* point here to start of buffer                    */ 
	  chunk_needed=here_position=here_chunkNumber=new_lock=0;  /* Reset all counters */
          lastfill_chunkNumber=-1;   /* fill buffer from the beginning                   */
      }
  }
      
  /* If we need more data,  and the last call did not return a new locked chunk */
  if( chunk_needed > k*npoint && new_lock == 0){
      int code;
      chunk_needed -= npoint;
      /* advance lastfill_chunkNumber (modulo kk)... */
      lastfill_chunkNumber = (lastfill_chunkNumber + 1) % kk;
      /* ...and fill the data in the chunk */
      code = fill_chunk(chunk + lastfill_chunkNumber,npoint);
      switch (code) {
      case 0:  /* no more data */
	  return 0;
      case 1:  /* entering new locked set */
	  new_lock = 1;
	  break;
      case 2:  /* continuing a locked set */
	  count_chunks++;
	  count_tinlock+=((double)npoint/SRATE);
	  break;
      default: /* unrecognized code */
	  GR_start_error("get_more_data()",rcsid,__FILE__,__LINE__);
	  GR_report_error("unrecognized code %d\n",code);
	  GR_end_error();
	  abort();
      }
      /* If data was put at beginning of buffer */
      if (lastfill_chunkNumber == 0) 
	  /* place a copy of it at the end to maintain continuity */
	  memcpy((void *)end_buffer,(void *)buffer,(size_t)npoint*sizeof(short));
  }

  /* compute the average spectrum for the kk chunks of buffered data */
  average_spectrum(chunk,kk,npoint,spec,new_lock);

  /* recalibrate (response only modified if needed) */
  recalibrate(response,npoint);

  /* the time of the returned data */
  *time = chunk[here_chunkNumber].time + (double)(here_position%npoint)/SRATE;

  /* copy the data to the returned array. */
  memcpy((void *)data,(void *)here,(size_t)npoint*sizeof(short));

  count_segments++;
  fprintf(fp_statistics,"Returning segment %d, based on chunks",count_segments);
  { int i;
  for (i=0;i<kk-1;i++)
      fprintf(fp_statistics," %d",chunk[i].counter);
  }
  fprintf(fp_statistics," %d\n",chunk[kk-1].counter);
  
  return out_code;
}


/* global variables for passing data from get_calibrated_data() */
extern double datastart;
extern float *n_inv_noise,*htilde,*pow_renorm,srate;
extern int npoint,new_lock,gauss_test;

int get_calibrated_data()
{
  void realft(float [], unsigned long, int);
#if(INSERT_CHIRP)
  void ins_chirp(int);
#endif
  
  static float *work,*spec,*response,*ave_spec;
  static short *data;
  static int first=1,npoint=NPOINT,num_sent=0;
  double fac=(double)npoint*SRATE;
  int i,cut,code,need=npoint-(POSTSAFETY+PRESAFETY);
  static int nomoredata=0;
  
  if (first) {
    first = 0;
    data = (short *)malloc(npoint*sizeof(short));
    work = (float *)malloc(npoint*sizeof(float));
    spec = (float *)malloc((npoint/2 + 1)*sizeof(float));
    ave_spec = (float *)malloc((npoint/2 + 1)*sizeof(float));
    response = (float *)malloc((npoint+2)*sizeof(float));
    if (data == NULL || work == NULL || spec == NULL
	|| ave_spec == NULL || response == NULL ) {
      GR_start_error("get_calibrated_data()",rcsid,__FILE__,__LINE__);
      GR_report_error("could not allocate memory\n");
      GR_end_error();
      abort();
    }
  }
  
  /* return 0 if the required number of segments (or more!) have been sent */
  /*  fprintf(stderr,"Just sent segment number %d\n",num_sent); */
  if ((num_sent++)>=DATA_SEGMENTS) return 0;
  if (nomoredata) return 0;
  
  code = get_more_data(need,data,ave_spec,response,&datastart);

  switch (code) {
  case 0:
      fprintf(fp_statistics,"Message from get_calibrated_data():\n");
      fprintf(fp_statistics,"  Last locked section had %f secs of data, ",MIN_INTO_LOCK*60.0 +
	      (float)(count_chunks)*npoint/srate);
      fprintf(fp_statistics," and %f secs discarded\n",discarded-lastdiscard);
      fprintf(fp_statistics,"No more data\n");
      fprintf(fp_statistics,"  Total locked data = %f secs, acquired lock %i times, ",
	      count_tinlock,count_locked);
      fprintf(fp_statistics," and analyzed %i segments\n",count_segments);
      fflush(fp_statistics); 
      nomoredata=1;
      return 0;
  case 1:
      new_lock = 1;
      break;
  case 2:
      new_lock = 0;
      break;
  default:
      GR_start_error("get_calibrated_data()",rcsid,__FILE__,__LINE__);
      GR_report_error("unrecognized code %d\n",code);
      GR_end_error();
      abort();
  }
  
  for (i=0;i<npoint;i++) work[i] = data[i];
  realft(work-1,npoint,1);
#if (RANDOMIZE)
  /* do everything BUT DC and Nyquist */
  for (i=1;i<npoint/2;i++) {
	int ir=i+i,ii=ir+1;
	double mag=sqrt(work[ir]*work[ir]+work[ii]*work[ii]);
	float phase=2.0*M_PI*ran2(&randomize);
	work[ir]=mag*cos(phase);
	work[ii]=mag*sin(phase);
  }
#endif
  product(htilde,work,response,npoint/2);
#if(INSERT_CHIRP)
  /* insert a chirrp if desired */
  ins_chirp(num_sent);
#endif
  
  /* lower cutoff frequency */
  cut = npoint*FLO/SRATE;
  if (cut<1) cut = 1;
  
  /* set n_inv_noise to zero at low frequencies and Nyquist */
  n_inv_noise[0] = 0;
  n_inv_noise[npoint/2] = 0;
  for (i=1;i<cut;i++) n_inv_noise[i] = 0;
  
  /* compute remaining n_inv_noise elements */
  for (i=cut;i<npoint/2;i++) {
    int ir=i+i,ii=ir+1;
    double re=response[ir],im=response[ii],tmp=ave_spec[i];
#if (NORM_CF)
    n_inv_noise[i] = 4.0/(fac*tmp*(re*re + im*im));
#else
    n_inv_noise[i] = 2.0/(fac*tmp*(re*re + im*im));
#endif
  }
  
# if (SPEC_TRUNC > 0) /* truncate time-domain version of n_inv_noise[]   */
  {
    float out0,norm = 2/(float)npoint;       /* normalization factor for iFFT */
    int spec_zero=(SPEC_TRUNC/2);            /* where to start zeroing        */
    for(i=0;i<npoint;i++) work[i]=0;         /* clear out work array          */
    for(i=0;i<npoint/2;i++) work[i+i]=sqrt(n_inv_noise[i]); /* fill the array */
    realft(work-1,npoint,-1);                /* iFFT it                       */
    for(i=spec_zero;i<npoint-spec_zero;i++) work[i]=0; /* truncate work[]     */
    realft(work-1,npoint,1);                 /* FFT to freq domain            */
    for(i=0;i<npoint/2;i++){
	out0=norm*work[i+i];               /* reconstruct n_inv_noise[]     */
	n_inv_noise[i]=out0*out0;
    }
    if ((cut = npoint*FLO/SRATE) < 1) cut = 1; /* low-frequency cutoff        */
    for (i=0;i<cut;i++) n_inv_noise[i]=0.0; /* clear low-frequency components */
    n_inv_noise[npoint/2] = 0;        /* make absolutely sure Nyquist is zero */
  }
# endif

#if (REMOVE_LINE_BINS)
  { 
    int harmonic;
    float freq;
    /* step through all frequency harmonics of 60 Hz. Note: freq=(srate*i)/npoint  */
    for (harmonic=1;harmonic<85;harmonic++) {
      float freqmin,freqmax;
      /* find the freq bins +- 0.5 N Hz above/below center */
      freqmin=harmonic*59.5;
      freqmax=harmonic*60.5;
      nullout(freqmin,freqmax,npoint,srate,n_inv_noise);
    }
    /* Now do the same for the violin resonances & other lines: */
    /* Bruce's list */
    nullout(  79.0,   81.0,npoint,srate,n_inv_noise);
    nullout( 109.0,  110.0,npoint,srate,n_inv_noise);
    nullout( 139.0,  140.0,npoint,srate,n_inv_noise);
    nullout( 245.0,  246.0,npoint,srate,n_inv_noise);
    nullout( 487.0,  490.0,npoint,srate,n_inv_noise);
    nullout( 499.0,  501.0,npoint,srate,n_inv_noise);
    nullout( 571.0,  572.0,npoint,srate,n_inv_noise);
    nullout( 576.0,  585.0,npoint,srate,n_inv_noise);
    nullout( 592.0,  602.0,npoint,srate,n_inv_noise);
    nullout( 603.0,  607.0,npoint,srate,n_inv_noise);
    nullout( 998.0, 1001.0,npoint,srate,n_inv_noise);
    nullout(1155.0, 1159.0,npoint,srate,n_inv_noise);
    nullout(1208.0, 1214.0,npoint,srate,n_inv_noise);
    nullout(1740.0, 1750.0,npoint,srate,n_inv_noise);
    nullout(3500.0, 3520.0,npoint,srate,n_inv_noise);
    
    /* Stan's list */
    nullout(  77.0,   83.0,npoint,srate,n_inv_noise);
    nullout( 105.0,  115.0,npoint,srate,n_inv_noise);
    nullout( 569.0,  607.0,npoint,srate,n_inv_noise);
    nullout(1138.0, 1214.0,npoint,srate,n_inv_noise);
    nullout(1707.0, 1821.0,npoint,srate,n_inv_noise);
    nullout(3500.0, 3520.0,npoint,srate,n_inv_noise);
  }
#endif /* (REMOVE_LINE_BINS) */
  
  /* compute outlier statistic */
  gauss_test = is_gaussian(data,npoint,-2048,2047,0);
  
  /* compute power statistic---don't include DC term */
  get_spectrum(data,spec,npoint);
  *pow_renorm = 0;
  for (i=1;i<npoint/2;i++) *pow_renorm += spec[i]/ave_spec[i];
  *pow_renorm *= 2.0/npoint;
  return 1;
}


#if(INSERT_CHIRP)
/* routine to determine if a chirp is to be inserted, and to insert it */
void ins_chirp(int segment)
{
  void realft(float *, unsigned long, int);
  static FILE *fpinsert,*fpinslog;
  static double instime=0;
  static float *chirp0,*chirp1,m1,m2,invMpc,c0,c1,phase;
  static int first=1,end=0,npoint=NPOINT;
  int offset;
  
  /* no more chirps to insert */
  if (end) return;
  if (first) {
    first = 0;
    
    /* open the insert.ascii file for input */
    fpinsert = grasp_open("GRASP_INSERT","insert.ascii","r");
    
    /* open the insert.log file for output */
    fpinslog = grasp_open("GRASP_INSERT","insert.log","w");
    
    /* allocate memory to chirp arrays */
    chirp0 = (float *)malloc(npoint*sizeof(float));
    chirp1 = (float *)malloc(npoint*sizeof(float));
  }
  
  /* scan through the file until the next chirp is found */
  /* note: assume that only one chirp will be present in any data segment! */
  while (instime<datastart) {
    
    float tc,scale=2*HSCALE*M_SOLAR/MPC;
    int i,code,chpts;
    
    /* read the next injected chirp time */
    code = fscanf(fpinsert,"%lf %f %f %f %f\n",
		  &instime,&m1,&m2,&invMpc,&phase);
    /* if we have reached the end of the file */
    if (code==EOF) {
      end = 1;
      fclose(fpinsert);
      fclose(fpinslog);
      free(chirp0);
      free(chirp1);
      return;
    }

    /* If injection is well before the segment start, try again */
    if (instime<(datastart-2*npoint/srate)) continue;
    
    /* coefficients of the injected chirp */
    c0 = cos(phase);
    c1 = sin(phase);
    
    /* construct the chirp to be injected */
#if (INJECT_TIME_REVERSE)
    make_retlifs(m1,m2,chirp0,chirp1,FLO,npoint,srate,&chpts,&tc,4000,4);
#else
    make_filters(m1,m2,chirp0,chirp1,FLO,npoint,srate,&chpts,&tc,4000,4);
#endif
    for (i=0;i<chpts;i++) {
      chirp0[i] *= scale;
      chirp1[i] *= scale;
    }
    for (i=chpts;i<npoint;i++) chirp0[i] = chirp1[i] = 0;
    realft(chirp0-1,npoint,1);
    realft(chirp1-1,npoint,1);
  }
  /* compute offset; return if after the end of the data segment */
  if (instime>(datastart+npoint/srate)) return;
  offset = srate*(instime - datastart);
  
  /* inject the chirp */
  freq_inject_chirp(c0,c1,offset,invMpc,chirp0,chirp1,htilde,npoint);
  
  /* write an entry into the log file */
  fprintf(fpinslog,"%d %d %f %f %f %f %f\n",
          segment,offset,instime,m1,m2,invMpc,phase);
  fflush(fpinslog);
  return;
}
#endif
