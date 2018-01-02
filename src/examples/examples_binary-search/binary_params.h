/* GRASP: Copyright 1997,1998  Bruce Allen */
/* Parameters that describe the size of the segments of data analyzed                     */
#define NPOINT    262144       /* The size of our segments of data (~26 seconds)          */
#define CHIRPLEN   24000       /* length of longest allowed chirp 23293 points            */
#define PRESAFETY  65536       /* ignore PRESAFETY points at beginning of correlation     */
#define POSTSAFETY (PRESAFETY + CHIRPLEN)  /* ignore POSTSAFETY at end of correlation     */
#define SPEC_TRUNC PRESAFETY
/* #define SPEC_TRUNC 0 */

/* Parameters that describe the data aquisition */
#define MIN_INTO_LOCK 3.0     /* Number of minutes to skip into each locked section       */
#define DATA_SEGMENTS 16384   /* largest number of data segments to process               */

/* Parameter that defines format of data being used                                       */
/* 0: old format data                                                                     */
/* 1: new (frame) format data                                                             */
/* 2: simulated data                                                                      */
#define DATA_FORMAT 1
#define RANDOMIZE 0           /* Randomizes the phases                                    */

/* Parameters that describe the interferometer */
#define FLO 120.0             /* The low frequency cutoff for filtering (Hz)              */
#define HSCALE 1.e21          /* A convenient scaling factor; results independent of it   */
#define SRATE 9868.4208984375 /* sample rate in HZ of IFO_DMRO channel                    */

/* Parameters that describe the filter bank and the filter signals */
#define NUM_TEMPLATES 687     /* number of templates read from an ascii file              */
#define STORE_TEMPLATES 0     /* 0: slaves recompute templates. 1: slaves save templates. */
#define NSIGNALS 8            /* number of signal values computed for each template       */
#define THRESHOLD 5.0         /* threshold above which a splitup veto test is done        */
#define PBINS 20              /* Number of frequency bands for the r^2 discriminator      */

/* Parameters controlling optional MPI and MPE code */
#define DEBUG_COMM 0          /* print lots of statements to help debug communication     */
#define DBX 0                 /* start up dbx windows for each process                    */
#define ISMPE 0               /* 0: no mpe calls.  1: mpe calls                           */
#define SMALL_MPELOG 1        /* 0: detailed MPE logging. 1: minimal MPE logging          */
#define RENICE 1              /* Priority. 1: processes run NICED  0: not NICED           */
#define NPERSLAVE 15           /* The number of segments analyzed in parallel by slave     */

/* Parameters for analysing a particular segment and a particular template */
#define PART_TEMPLATE -1      /* template to be investigated (numbered 0,1,...)           */
#define COMPARE 0             /* 0: records impulse offset. 1: records peak offset        */

/* parameters for type of filtering */
#define INSERT_CHIRP 0        /* 0: don't insert a chirp.  1: insert chirp into data      */
#define REMOVE_LINE_BINS 0    /* implement Whitcomb's idea of dropping bins near lines    */
#define INJECT_TIME_REVERSE 0 /* inject time-reversed chirps                              */
#define REVERSE_FILTERS 0     /* analyze data with time-reversed filters                  */
#define WINDOW 2              /* power spectrum: 0=rectangular 1=Hann 2=Welch 3=Bartlett  */
#define NORM_CF 1             /* 0: GRASP snr/distance norm  1: Cutler/Flanagan norms     */
#define TWO_PHASE_R2 1        /* 1: Use two-phase r^2 test, 0: use single-phase test      */

/* Utility parameters  */
#define HEADER_COMMENT_SIZE 16384  /* bytes in the comment part of signal.header file     */
