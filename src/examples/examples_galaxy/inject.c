/* GRASP: Copyright 1997,1998  Bruce Allen */
/* Program inject.c

   Produces a list of random Galactic NS-NS binary inspiral parameters for
   injection into the interferometer data to test the data analysis software.

   The inspiral waveforms correspond to two NS companions, each with a mass
   distribution that is uniform between two cutoffs MLO and MHI
   [for example see L S Finn, Phys Rev Lett 73 1878 (1994)].

   The amplitude distribution corresponds to the (time dependent) model
   described for the GRASP routine mc_chirp().  The initial phase is uniformly
   distributed.  The variables LAT, LON, and ARM correspond to the latitude,
   longitude, and arm orientation of the detector, and are required for the
   model.

   The injection time is either at a fixed intervals INV_RATE (FIXED=1),
   or at random intervals corresponding to a Poisson process with
   inverse rate INV_RATE (FIXED=0).
   Injection times are between the start and the end of the data run specified
   by the environment variable GRASP_DATAPATH.  The start and end times of this
   data run are obtained from code resembling that in program locklist.c.
   If two chirps potentially occur within the same data segment (of length
   NPOINT points), a warning message is printed.

   The results are output to stdout in a list containing the arrival time
   (double), the two masses (floats), the amplitude---inverse Mpc distance
   (float), and the initial phase (float), separated by spaces.  This is the
   same format as required for the file insert.ascii which is read by the
   binary_get_data() routine in the binary search code.
*/
#include "grasp.h"

#define OFFSET 15.0   /* offset in secs of injected chirps arrival times */
#define SEED -101     /* initial seed value for random #, <0 */
#define LAT  34.1667  /* detector latitude in degrees North */
#define LON 118.133   /* detector longitude in degrees West */
#define ARM 180.0     /* detector arm orientation in degrees CCW from North */
#define MLO 1.29      /* low NS mass limit in solar masses */
#define MHI 1.45      /* high NS mass limit in solar masses */
#define INV_RATE 30.0 /* inverse of event rate in seconds */
#define FIXED 1       /* 0: Events occur in Poisson-distributed intervals
			 1: Events occur at the fixed rate RATE */
#define NPOINT 262144 /* number of points in each data segment (used to warn
			 when two events may exist in the same segment) */
#define KNOWN_START_END 1 /* use the starting times below rather than from data */
#define TSTART 784880277
#define TEND   785388428
float ran1(long *);

int known_time(double *time,long *seed){
  static int first=1;
  if (first) {
    *time=TSTART+OFFSET;
    first=0;
  }
  else {
    if (FIXED)
      *time+=INV_RATE;
    else
      *time+=INV_RATE*log(ran1(seed));
  }
  return (*time<TEND);
}


int main()
{
  float local_sidereal_time(time_t, float);
  void mc_chirp(float, float, float, long *, float *, float *, float *);
  int next_time(double *, long *);
  double time;
  long seed=SEED;
  
#if (KNOWN_START_END)
  while (known_time(&time,&seed))
#else
  while (next_time(&time,&seed))
#endif
    { 
      float lst,m1,m2,c0,c1,phase,invMpc;
      
      /* compute the local sidereal time in seconds */
      lst = 3600*local_sidereal_time((time_t)time,LON);
      
      /* obtain the random chirp parameters */
      mc_chirp(lst,LAT,ARM,&seed,&invMpc,&c0,&c1);
      
      /* the random phase in radians */
      phase = atan2(c1,c0);
      
      /* the mass of the first and second NSs */
      m1 = MLO + (MHI - MLO)*ran1(&seed);
      m2 = MLO + (MHI - MLO)*ran1(&seed);
      
      /* print the parameters */
      printf("%f %f %f %f %f\n",time,m1,m2,invMpc,phase); 
    }  
  return 0;
}


/*
  Calculates the next time of a binary inspiral and returns 1 if this time
  occurs before the end of the data run or 0 if it does not.
  */
int next_time(double *time, long *seed)
{
  float ran1(long *);
  
  static int first=1;
  static double end;
  static float srate;
  float dt;
  
  if (first) { /* obtain the start and the end of the data run on first call */
    
    time_t begin;
    float st,et,slowrate;
    int soff,sblk,eoff,eblk;
    FILE *fplock;
    
    first = 0;
    
    /* read ld_mainheader to get begin time */
    {
      
      struct ld_mainheader mh;
      struct ld_binheader bh;
      short *datas=NULL;
      FILE *fpifo;
      int alloc=0,n=0;
      fpifo = grasp_open("GRASP_DATAPATH","channel.0","r");
      read_block(fpifo,&datas,&n,&st,&srate,1,&alloc,1,&bh,&mh);
      fclose(fpifo);
      begin = (time_t)mh.epoch_time_sec;
    }      
    /* start of the first locked segment */
    fplock = grasp_open("GRASP_DATAPATH","channel.10","r");
    find_locked(fplock,&soff,&sblk,&eoff,&eblk,&st,&et,&slowrate);
    *time = (double)(begin) + (double)(soff/srate) + OFFSET;
    
    /* end of the last locked segment */
    while (find_locked(fplock,&soff,&sblk,&eoff,&eblk,&st,&et,&slowrate));
    fclose(fplock);
    end = (double)(begin) + (double)(et + eoff/srate);
    
    /* print the start and end times and the duration of the run */
    fprintf(stderr,"start: %f, end: %f, duration: %f\n",*time,end,end - *time);
    
  }
  if (FIXED) /* fixed rate intervals */
    dt = INV_RATE;
  else /* intervals for a Poisson process */
    dt = -INV_RATE*log(ran1(seed));
  
  if (dt*srate<NPOINT) /* print warning message when chirps too close */
    fprintf(stderr,"Warning: potentially two chirps in same segment\n");
  
  /* increment the time */
  *time += dt;
  
  if (*time<end) return 1;
  else return 0;
}
