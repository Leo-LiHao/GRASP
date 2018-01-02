/*GRASP: Copyright, 1997,1998 Bruce Allen */

#define DATADIM 65536    /* size of segment  (must be power of 2)*/
#define NDIM   4096      /* size of the subsegment to be transformed (must be power of 2)*/
#define PDIM 512         /* dimension of time-freq map; (max value NDIM/4)  */
#define PRESAFETY DATADIM/8  /* number of points to ignore at beginning of correlation   */
#define POSTSAFETY DATADIM/8 /* number of points to ignore at end of correlation        */
#define FDIM (NDIM/(4*PDIM))  /* number of bins to average over in frequency space*/
#define TDIM (NDIM/PDIM)   /* number of bins to average over in time space */

/* THE TYPE OF NOISE */
#define NOISE_WHITE 1    /* white noise */
#define NOISE_LIGO_INI 2 /* initial ligo noise curve */

/* THE TYPE OF SIGNAL */
#define INSERT_INSPIRAL 1     /* the inspiral waveform */
#define INSERT_QUAS_PER 2     /* a waveform with power law increase of freq and amp*/
#define INSERT_COALESCENCE 3 /* the full coalescence waveform */

#define RESCALE_FACTOR 2.0 /* used to decide the rescale number*/
#define DEBUG  0          /* For debugging purposes */
#define DEBUG1 0          /* temprorary debugging */
#define DEBUG2 0           /* temprorary debugging */
#define NOISE_TYPE NOISE_LIGO_INI  /* To determine the noise_type */

#ifndef SRATE
#define SRATE 9868.4208984375 /* sample rate in HZ of IFO_DMRO channel                    */
#endif

/* Function declarations */
void fill_data_with_signal(int, float *, float *, float);
void get_coalescence(float *, int ,  float , float ,int *);
int get_time_series_data();
int gethtilde();
void getshf(int,float *,float);
void gettfparameters();
void master();
float mygasdev(long *);
void noise_gau_col(long *,unsigned long,float *, float *);
void noise_gau_col_fr(long *,unsigned long , float *, float *);
void normalizehtilde(float *, int,float * );
void over_whiten_filter(float *, long , float *);
void picturerawprint(float **);
void realft(float *, unsigned long, int);
void slave();
void timstat(int, FILE *, float *htilde);
void whiten_filter(float *, long , float *);

/* Structure definitions */
typedef struct 
{
    int signaltype;        /* the type of signal to be used */
    float m1;              /* the mass of one of the stars; used if signaltype = 1 */
    float m2;              /* the mass of the other star; used if signaltype = 1 */  
    float pind;            /* the exponent for the power law freq. increase;used if signaltype = 2 */  
    float ampind;          /* the exponent for the power law amplitude. increase;used if signaltype=2 */  
    float timfrac;         /* the fraction of time for which the signal lasts in each subsegment
                              ;used if signaltype=2 */
    float freqfrac;        /* the bandwidth of the signal as a fraction of the sampling rate;
                              used if signaltype=2 */
    float cons_noise_pow;  /* an arbitrary scaling factor for the power spectrum */
    float snr;             /* the signal to noise ratio */
    int signaloffset;      /* the offset at which the signal begins in each subsegment */
    int addsignal;         /* flag; TO ADD signal 1; if signal is not to be added 0 */
} struct_signalparameters;
