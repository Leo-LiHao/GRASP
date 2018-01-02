/* GRASP: Copyright 1997,1998  Bruce Allen */

#ifndef _GRASP_H
#define _GRASP_H

/*  We Require that all code be Posix-compliant */
#ifndef _POSIX_SOURCE /* Avoid redefinition */
#define _POSIX_SOURCE
#endif

#ifdef __cplusplus /* Explicit C linkage is required in C++ code */
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include <stdarg.h>
    
#ifndef M_PI
#define M_PI            3.14159265358979323846
#endif
    
#ifndef M_PI_2
#define M_PI_2          1.57079632679489661923
#endif

/* Astrophysical constants from Particle Physics Data Book              */
#define C_LIGHT   (2.99792458e10)        /* speed of light [cm/sec]     */
/* quantity in Particle Physics Data Book, divided by two               */
#define M_SOLAR   (1.47662504e5)         /* solar mass [cm]             */
#define TSOLAR    (M_SOLAR/C_LIGHT)      /* solar mass [seconds]        */
#define MPC	  (3.0856775807e24)      /* one mega parsec [cm]        */
#define HUBBLE    (3.227296e-18)         /* sec^-1 h_100                */

/* difference between UTC and GPS time 315964811 = 
3600 sec/hour x 24 hours/day x (365 days/year x 8 years + 366 days/year x 2 years+ 5 days) + 11 leap seconds */
#define UTCTOGPS  315964811              
#define MAX_PHS_TERMS   15

/* Structure used in the November 1994 40-meter data format */
struct ld_binheader {
  float elapsed_time;
  float datarate; };

/* Structure used in the November 1994 40-meter data format */
/* Note that we have changed long->int from the original definition
   because we need 4-byte quantities.  On some machines (dec alpha)
   they came across as 8-bytes
*/
struct ld_mainheader{
  int chunksize;
  int filetype;
  int epoch_time_sec;
  int epoch_time_msec;
  int tod_second;   /* seconds after minute, 0-61 for leap second */
  int tod_minute;   /* minutes after hour 0-59 */
  int tod_hour;     /* Hour since midnight 0-23 */
  int date_day;     /* 1-31 */
  int date_month;   /* 0-11 is Jan- Dec */
  int date_year;    /* years since 1900 */
  int date_dow;     /* days since Sunday 0-6 */
  int sub_hdr_flag; };

/* average optical arm-cavity length (meters) for 40-meter in November 1994 */
#define ARMLENGTH_1994 38.25

struct Template {
	int num;	/* identifies which template */
	float f_lo;	/* start (low) freq of template, sec^-1 */
	float f_hi;	/* end freq of template, sec^-1 */
	float tau0;	/* Newtonian time to coalescence, sec */
	float tau1;	/* First post-Newtonian correction to tau0 */
	float tau15;	/* 3/2 PN correction */
	float tau20;	/* second order  PN correction */
	float pha0;	/* Newtonian phase to coalescence, radians */
	float pha1;	/* First post-Newtonian correction to pha0 */
	float pha15;	/* 3/2 PN correction */
	float pha20;	/* second order  PN correction */
	float mtotal;	/* total mass, in solar masses */
	float mchirp;	/* chirp mass, in solar masses */
	float mred;	/* reduced mass, in solar masses */
	float eta;	/* reduced mass/total mass, dimensionless */
	float m1;	/* smaller of the two masses */
	float m2;	/* larger of the two masses  */
};

/* variables to describe the scope of the search (the family of templates) */
struct  Scope {
	int n_tmplt;	/* computed number of templates required for search */
	float m_mx;	/* maximum mass in our search */
	float m_mn; 	/* minimum mass in our search */
	float theta; 	/* angle (rad) orientation in template space */
	float dp; 	/* spacing between templates  */
	float dq; 	/* spacing between templates  */
	float f_start; 	/* frequency at the time of arrival  */
	struct Template *templates; /* pointer to list of templates */
};


struct Raw_data_header {
	int n;		/* number of (valid) data points */
	int block;	/* which block of data from large set */
	float t_sample;	/* sample time, in sec */
	float f_sample;	/* sample freq, in Hz */
	double t_start;	/* time stamp of the first data point */
	short *begin;	/* location in memory of the first data point */
};


struct Swept_sine_header {
	int n;		/* number of points */
	int block_lo;	/* valid from this data block... */
	int block_hi;	/* to this one (inclusive) */
	float *freq;	/* pointer to frequencies (Hz) */
	float *real;	/* pointer to real voltage response */
	float *imag;	/* pointer to imaginary voltage response */
};	


/* Calibrated data always in freq domain! */
struct Cal_data_header {
	int n;		/* number of (valid) data points */
	int block;	/* which block of data from large set */
	float t_sample;	/* sample time, in sec */
	float f_sample;	/* sample freq, in Hz */
	double t_start;	/* time stamp of the first data point */
	float f_lo;	/* start frequency */
	float f_hi;	/* end frequency */
	int format;	/* how data packed (see above) */
	float *begin;	/* location in memory of the first data point */
};


struct Signal {
	float t_arr;	/* arrival time maximizing snr, sec */
	float pha_arr;	/* arrival phase maximizing snr, radians */
	float snr_max;	/* maximum signal to noise ration */
};

/* used in multitaper for removing spectral lines */
struct removed_lines{
        int index;      /* bin of spectral line removed */
	float fvalue;   /* value of F-test for corresponding line */
	float re;       /* real part of Fourier coefficient C removed */
	float im;       /* imaginary part of Fourier coefficient C removed */
	float dcdbr;	/* real part of dC/dbin */
	float dcdbi;	/* imag part of dC/dbin */
	float d2cdb2r;	/* real part of d^2C/dbin^2 */
	float d2cdb2i;	/* imag part of d^2C/dbin^2 */
};


int get_data(FILE*,FILE*,float*,int,short*,int*,float*,int);
int get_data2(FILE*,FILE*,float*,int,short*,int*,float*,int);
void normalize_gw(FILE*,int,float,float*);
void template_grid(struct Scope *Grid);
void plot_template(char*,struct Scope,int,int);
void calibrate(FILE * ,int num,float *complex,float srate,int method,int order);
int find_locked(FILE*,int*,int*,int*,int*,float*,float*,float*);
int find_locked2(FILE*,int*,int*,int*,int*,float*,float*,float*);
float template_area(struct Scope*);
void read_sweptsine(FILE*,int*,float**,float**,float**);
int read_block(FILE*,short**,int*,float*,float*,int,int*,int,struct ld_binheader*,struct ld_mainheader*);
void tau_of_mass(double,double,double,double*,double*);
int m_and_eta(double tau0,double tau1,double* mtot,double* eta,double mmin,double mmax,double pifreq);
void avg_spec(float *data,float *average,int npoint, int *reset,float srate,float decaytime,int windowtype,int overlap);

void ratio(float *c,float *a, float *b,int ncomplex);
void product(float *c,float *a, float *b,int ncomplex);
void productc(float *c,float *a, float *b,int ncomplex);
void ratio(float *c,float *a, float *b,int ncomplex);
void correlate(float *c,float *a,float *b,float *s,int n);
int phase_frequency(float,float,float,float,int,float*,float,float,float*,
                    float,float**,float**,int*,int*,float*,int);
int chirp_filters(float,float,float,float,int,float*,float,float,float*,
                    float,float**,float**,int*,int*,float*,int);
void make_filters(float m1,float m2,float *ch1,float* ch2,float
	fstart,int steps_alloc,float srate,int *filled,float *t_coal,int,int order);
void binner(float *input,int ninput,double *bins,int nbins,
            double firstbin,double width,int *under,int *over);
void binshort(short *input,int ninput,double *bins,int offset);
void clear(float*,int,int);
void graph(float *array,int n,int spacing);
void graph_double(double *array,int n,int spacing);
void graph_short(short *array,int n);
void splitup(float *c,float *a,float *s,int n,float total,int segments,int *indices);
void sgraph(short *array,int n,char* ext,int filenumber);
void audio(short *array,int n);
void sound(short *array,int n,char *,int filenumber);
int is_gaussian(short *array,int n,int min,int max,int print);
void freq_inject_chirp(float c0,float c90,int offset,float invMpc,float* ch0tilde,float* ch1tilde,float* htilde,int n);
void time_inject_chirp(float c0,float c90,int offset,float invMpc,float *chirp0,float *chirp1,float* data,float* response,float* work,int n);
float splitup_freq(float c0,float c1,float *chirp0, float *chirp1,float norm,
		float* twice_inv_noise,int npoint,int offset,int nbands,int* indices,float* stats,
			float* working,float* htilde);
float splitup_freq2(float c0,float c1,float *chirp0, float *chirp1,float
		norm, float* twice_inv_noise,int npoint,int offset,int
		nbands,int* indices,float* stats,
		float* working,float* htilde);
float splitup_freq3(float c0,float c1,float *chirp0, float *chirp1,
		    float norm, float* twice_inv_noise,int npoint,
		    int offset,int nbands,int* indices,float* stats, 
		    float* working,float* htilde); 
float splitup_freq4(float c0,float c1,float *chirp0, float *chirp1,
		    float norm, float* twice_inv_noise,int npoint,
		    int offset,int nbands,int* indices,float* stats, 
		    float* working,float* htilde); 
float splitup_freq5(float c0,float c1,float *chirp0, float *chirp1,
		    float norm, float* twice_inv_noise,int npoint,
		    int offset,int nbands,int* indices,float* stats, 
		    float* working,float* htilde);
void avg_inv_spec(float flo,float srate,int npoint,double decay,double *norm,
                  float *htilde, float* mean_pow_spec,float* twice_inv_noise);
void orthonormalize(float* ch0tilde,float* ch1tilde,float* twice_inv_noise,int npoint,float* n0,float* n1);

void find_chirp(float* htilde,float* ch0tilde,float* ch1tilde,float* twice_inv_noise,float n0,float n1,
		float* output0,float* output1,int npoint,int chirplen,int* offset,float* snr_max,
		float* lin0,float* lin1,float *var);

void sp_filters(float m1, float m2, float *ch0, float *ch90, float fstart,
    int n, float srate, float f_c, float t_c, int order);

FILE * grasp_open(const char *environment_variable,const char *shortpath,const char *mode);

void detector_site(char *,int,float [],char *,char *,char *);
void noise_power(char *,int,float,double *);
void whiten(char *,int,float,double *);
void overlap(float *,float *,int,float,double *);
void normalize(float *);
void monte_carlo(int,int,int,int,float,float,float,float,double *,double *,
                 double *,double *,double *,float *,float *,int *);
void get_IFO12(FILE *,FILE *,FILE *,FILE *,int,float *,float *,float *,float *);
void simulate_sb(int,float,float,float,float,double *,double *,double *,
		   float *,float *,int *);
void simulate_noise(int,float,double *,double *,float *,int *);
void combine_data(int,int,float *,float *,float *);
int test_data12(int n,float *,float *);
void analyze(int,float *,float *,int,float,float,float,float,double *,
	     double *,double *,int,int,double *,double *,double *,double *);
void extract_noise(int,int,float *,int,float,double *,double *);
void extract_signal(int,float *,float *,int,float,double *,double *,double *);
void optimal_filter(int,float,float,float,double *,double *,double *,double *);
double calculate_var(int n,float,float,float,float,float,double *,double *,
                     double *);
void prelim_stats(float,float,double,double);
void statistics(float *,int,int);

void remove_spectral_lines(float *data,int npoints,int padded_length,
   float nwdt,int nwin, int max_lines, int maxpass, int *num_removed, 
   struct removed_lines *line_list,float *mtap_spec_init,float *mtap_spec_final,int dospec,int fimin, int fimax);

void slepian_tapers(int n, int nwin, double *el, float nwdt, double *tapers, double *tapsum);

void multitaper_spectrum(float *data, int npoints, int kind,
	int nwin, float nwdt, int inorm, float dt,
	float *ospec, float *dof, float *Fvalues, int klen,float *cest,int dospec);

void multitaper_cross_spectrum(float *data1,float *data2, int npoints, 
 			       int npadded, float dt, int nwin, float nwdt,
 			       double *spec12);

int fvalue_cmp(const void *l1, const void *l2);

int index_cmp(const void *l1, const void *l2);

int fget_channel(char *channelname, char *(*files)(),double *tstart,int npoint,short *location,double *srate);

void GRcalibrate(float *fri,int frinum,int num,float *complex,float srate,int method,int order);

void GRnormalize(	float *fri,	/* array to get frequency, real, imaginary parts.        */
			int frinum,	/* number of points in previous array (divisible by 3)   */
			int npoint,	/* number of points in time domain sample (power of 2)  */
			float srate,	/* sample rate in Hz of the gw signal stream            */
			float *response /* on return, contains complex response function        */
);

char *framefiles();

/* this is the input structure for fgetintput */
struct fgetinput {
	int nchan;					/* number of channels to get */
	char **chnames;					/* list of channel names to get */
	int npoint;					/* number of points to get from first channel listed */
	short **locations;				/* list of locations to put channel data into */
	char *(*files)();				/* function to call to get next non-empty file */
	int (*filedes)();				/* function to call to get next file descriptor */
	int inlock;					/* set to zero to get any contiguous data, non-zero for locked only */
	int seek;					/* set to zero to get data, set non-zero to seek forward */
	int calibrate;					/* do we want calibration f,r,i? */
        char *datatype;                                /* character strings indicating data type in each channel */
};

/* this is the output structure for fgetintput */
struct fgetoutput {
	double tstart;					/* time of first point in output */
	double tstart_gps;				/* time of first point in output (GPS time)*/
	double srate;					/* sample rate of first requested channel */
	int *npoint;					/* number of points returned in each channel (array) */
	int *ratios;					/* ratios of sample rates channel[0]/channel[i]      */
	int discarded;					/* number of points discarded (out of lock, frames not continuous) */
	double tfirst;					/* time of first data point in the run (GPS time) */
	double tfirst_gps;				/* time of first data point in the run */
	double dt;					/* tstart-tfirst */
	double lostlock;				/* last loss of lock after previous */
	double lostlock_gps;				/* last loss of lock after previous */
	double lastlock;				/* last regain of lock  */
	double lastlock_gps;				/* last regain of lock  */
	int returnval;					/* integer value returned */
	float *fri;					/* pointer to freq, real, imag part of calibration info */
	int frinum;					/* number of frequency values in calibration info */
	int tcalibrate;					/* time of last calibration (seconds) */
	int tcalibrate_gps;				/* time of last calibration (seconds) */
	int locklow;					/* min value at which we are still in lock */
	int lockhi;					/* max value at which we are still in lock */
	char *filename;					/* name of most recently opened frame file (or NULL) */
        char *slow_names;                               /* names of the slow channels packed into one fast channel SLOW */
};

int fget_ch(struct fgetoutput *fgetoutput,struct fgetinput *fgetinput);
struct qnTemplate {
	int num;
	float freq;
	float qual;
};

struct qnScope {
	int n_tmplt;
	float freq_min;
	float freq_max;
	float qual_min;
	float qual_max;
	struct qnTemplate *templates;
};

void qn_eigenvalues(float eigenvalues[], float a, int s, int l, int m);
void sw_spheroid(float *re, float *im, float mu, int reset,
		 float a, int s, int l, int m, float eigenvalues[]);
int qn_ring(float iota, float beta,
            float eps, float M, float a, int l, int m,
            float du, float atten, int max,
	    float **plusPtr, float **crossPtr);
int qn_qring(float psi0, float eps, float M, float a,
             float du, float atten, int max, float **strainPtr);
int qn_filter(float freq, float qual,
              float du, float atten, int max, float **filterPtr);
void qn_normalize(float *u, float *q, float *r, int n, float *norm);
void find_ring(float *h, float *u, float *r, float *o,
               int n, int len, int safe, int *off,
               float *snr, float *mean, float *var);
void qn_inject(float *strain, float *signal, float *response, float *work,
               float invMpc, int off, int n, int len);
void qn_template_grid(float dl, struct qnScope *grid);

/* types for all GRASP error handling functions */
#ifndef __CINT__ /* for ROOT */

	typedef void (*GR_start_error_type)(const char* function, const char* rcsid,
	                                    const char* file,     const long line);
	typedef void (*GR_report_error_type)(const char* format, va_list arglist);
	typedef void (*GR_end_error_type)(void);

/* GRASP's interface to the routines pointed to by the above pointers */

	void GR_start_error(const char* function, const char* rcsid,
	                    const char* file,     const long line);
	void GR_report_error(const char* format, ...);
	void GR_end_error(void);

	int  GR_is_reporting(void);

	void GR_set_error_handlers
		(
		GR_start_error_type  error_start,
	    GR_report_error_type error_report,
	    GR_end_error_type    error_end
	    );
	void GR_get_error_handlers
		(
		GR_start_error_type*  error_start,
	    GR_report_error_type* error_report,
	    GR_end_error_type*    error_end
	    );

	void GR_set_errors_enabled(int on_off);
	int  GR_errors_enabled(void);

/* The GRASP default error package */

	void GR_default_start_error(const char* function, const char* rcsid,
	                            const char* file,     const long line);
	void GR_default_report_error(const char* format,   va_list arglist);
	void GR_default_end_error(void);

	void GR_restore_default_handlers(void);

	void  GR_set_error_file(FILE* new_error_file);
	FILE* GR_get_error_file(void);

	int         GR_set_error_file_name(const char* new_error_file_name,
	                                   const int   erasefile);
	const char* GR_get_error_file_name(void);

#endif

/* The prototypes for template match fitting routines */

void bisect(float [],float [],float [],float,int *,float [],
	    float [],float,float,float,float,float,float [],
	    float,float,float,float [],int,float [],float[],int);
int parab(float [],float [],float [],int,float []);
int cubic(float [],float [],float [],int,float []);
float compute_match(float,float,float [],float [],float,float [],
		float,float,float,int,float,int,int);
int match_cubic(float,float,float,int,float,float,float,char *,
		float *,float *,float *,float []);
int match_parab(float,float,float,int,float,float,float,char *,
		float *,float *,float *,float []);


/*
 *
 * GENERAL UTILITIES
 *
 */


/* computes the local sidereal time (decimal hours) at a given time
   at a given longitude (degrees West) */
float local_sidereal_time(time_t time, float longitude);


/*
 *
 * MONTE CARLO ROUTINES
 *
 */


/* produces a random (effective) amplitude and phase coefficients of a
   Galactic NS-NS binary at a given local sidereal time (seconds) for
   a detector at a given latitude (radians N) and with a given arm
   orientation (radians CCW from N) */
void mc_chirp(float time, float latitude, float orientation, long *seed,
	      float *invMpc, float *c0, float *c1);

/* computes the beam pattern functions, plus and cross, for some specified
   angles, theta, phi, and psi (in radians) */
void beam_pattern(float theta, float phi, float psi,
		  float *plus, float *cross);

/* generates a random sky position of a NS-NS binary */
void sky_position(float u1, float u2, float u3,
		  float *rad, float *alpha, float *delta);

/* converts from equatorial to horizon coordinates for a given
   time (sidereal seconds) and latitude (rad)
   note: alpha, delta, azi, and alt all in radians */
void equatorial_to_horizon(float alpha, float delta, float time, float lat,
			   float *azi, float *alt);

/* converts from galactic coordinates to equatorial coordinates
   note: l, b, alpha, and delta all in radians */
void galactic_to_equatorial(float l, float b, float *alpha, float *delta);

/* Prototypes for template placement routines */
 struct cubic_grid {
   int n;
   float m_mn;
   float m_mx;
   float dm;
   float match;
   float angle;
   int order;
   float srate;
   float flo;
   float ftau;
   int detector;
   float ***coef;
 };
void generate_cubic(struct cubic_grid grid, char *detectors_file,
		    const char *outfile, const char *logfile);
int regenerate_cubic(char *detectors_file, const char *infile,
		     const char *outfile, const char *logfile);
int read_cubic(struct cubic_grid *grid, const char *infile);
void write_cubic(struct cubic_grid grid, const char *outfile);
int get_cubic(float m1, float m2, struct cubic_grid grid,
	      float *coef);
void free_cubic(struct cubic_grid grid);
void transform_cubic(struct cubic_grid *grid, float angle,
		     float match);
 struct tile {
   int flag;
   double x;
   double y;
   double dx;
   double dy;
   double r1;
   double r2;
   double theta;
   struct tile *next;
 };

 int tiling_2d(double *x_bound, double *y_bound, int npts,
	      int (*metric)(double , double , double *),
	      struct tile **tail, int *n_tiles);
int column(double *x_bound, double *y_bound, int npts, double left,
	   double width, int (*metric)(double , double , double *),
	   struct tile **tail, int *n_tiles);
int corner(double *x_bound, double *y_bound, int npts, int xsign,
	   int ysign, int (*metric)(double , double , double *),
	   struct tile *here, struct tile **tail, int *n_tiles);
int test_cycle(double *cycle, int npts, double *low, double *high);
int range(double *x_bound, double *y_bound, int npts, double x,
	  double *y_low, double *y_high);
int get_height(int (*metric)(double , double , double *), double x,
	       double y, double dx, double *dy);
int get_width(int (*metric)(double , double , double *), double x,
	      double y, double *dx);
int get_ellipse(int (*metric)(double , double , double *), double x,
		double y, double *r1, double *r2, double *theta);
int get_plot_bounds(double *x_bound, double *y_bound, int npts,
		    double angle, double magnification,
		    double internal_margin, double *x_corner,
		    double *y_corner, int *horizontal_pages,
		    int *vertical_pages);
int plot_list(double *x_bound, double *y_bound, int npts,
	      struct tile *head, int n_tiles, double angle,
	      double magnification, int plot_boundary,
	      int plot_tiles, int plot_ellipses, int plot_flags,
	      const char *psfile);
 struct chirp_space {
   float m_mn;
   float m_mx;
   float ftau;
   float angle;
   float match;
   int n_bound;
   double *x_bound;
   double *y_bound;
   struct cubic_grid grid;
   int n_templates;
   struct chirp_template *templates;
 };
 struct chirp_template {
   int flag;
   int num;
   double x;
   double y;
   double dx;
   double dy;
   double semimajor;
   double semiminor;
   double theta;
   double tau0;
   double tau1;
   double mtotal;
   double mchirp;
   double mred;
   double eta;
   double m1;
   double m2;
 };
void set_chirp_space(struct chirp_space space);
int chirp_metric(double x, double y, double *m);
void get_chirp_boundary(struct chirp_space *space);
int get_chirp_grid(struct chirp_space *space, const char *gridfile);
int get_chirp_templates(struct chirp_space *space);
int plot_chirp_templates(struct chirp_space space,
			 double magnification, int plot_boundary,
			 int plot_patches, int plot_ellipses,
			 int plot_flags, const char *psfile);

 /*                     
  *                     
  *                     TRANSIENT SOURCES ROUTINES
  *
  */
/* 
 * The differential equation for frequency for Lai-Shapiro core hang-up.
 */
void LS_freq_deriv(float t, float u[], float dudt[]);

/* 
 * integrates phase and frequency evolution equations for Lai-Shapiro core
 * hang-up.
 */
void LS_phas_and_freq(double Phi[], float u[], float fmax, float dt, 
                        int n_samples);
 
 /* structure to store physical parameters of Lai-Shapiro core hang-up */
 struct LS_physical_constants {
   float mass; 
   float radius; 
   float distance;
   float fmax;
   float inclination; 
   float Phi_0;
 };

 /*
  * generates detector stress waveform from Lai-Shapiro core hang-up.
  */
 void LS_waveform(float h[], struct LS_physical_constants phys_const,
            float sky_theta, float sky_phi, float polarization, float dt, 
          int n_samples);


/*
 *
 * TESTMASS ROUTINES
 *
 */

/* includes for the testmass package  S. Droz 1998 */
   
extern const int   kNumberOfFloats;       /* Number of floats allocated in one chunk */
enum testmass_errors  { kBhptNoError,     /* No error */
 					kBhptCantOpenFile,    /* Error opening file */
                    kBhptOutOfMemory,     /* Not enough memory to qperform operation */
                    kBhptUnknownMemory,   /* Don't know how much memory was allocated */
                    kBhptNotEnoughPoints, /* Not enough data points available */
                    kBhptCorruptFile,     /* Data file is corrupt */
                    kBhptStepTooSmall,    /* The stepsize is too small */
                    kBhptTooManySteps,    /* Too many steps needed */    
		    		kBhptNoDataRead,      /* You have not read any data in. */
		    		kBhptNoPhase,         /* The phase has not been calculated */
                    kBhptNoTime,          /* t(v) has not been calculated */
		    		kBhptFOutOfRange	  /* Illegal frequency requested */  };

/* Prototypes */
int integrateODE(float ystart[],  int nvar, float *x1, float x2, float eps, 
               float h1, float hmin, void (*derivs)(float, float [], float []));
int integrate_function(float vl, float vr, float vo, float (*f)(float ),
                        float **F, int number_of_points);
			
int read_real_data_file(const char *filename, float **x, float **y, int *number_of_points,
     int ReadX);
int read_modes(const  char *filename, float **x, float **ReA, float **ImA, 
               int *number_of_points, int *MaxL, int ReadX) ;     
float minustwoSlm(float theta, int l, int m);
		
int calculate_testmass_phase(float fo, float M, float **Phi);
int ReadData(char *filenameP,  char *filenameAlm, float **v, int *number_of_points);		 
void Clean_Up_Memory( float *Phase );
void Set_Up_Data( float *v, float *P, float *T, float *ReA, float *ImA,
                  int num_of_datapoints  );

int testmass_chirp(float m1, float m2, float theta, float phi,
               float *Phase, float f_start, float f_end, float *f_started,
               float *f_ended,  float dt, float **hplus, float **hcross, 
	       float **fre, int *number_of_points, int MaxL, int *modes);

float Get_Duration(float f1, float f2, float m1, float m2);
float Get_Fmax(float m1,float m2);


void merger_dist(double *deff, double *z, double *Vc, double m1_z,
		 double m2_z, double snr, double S_h[], int npoint,
		 double srate,double h100);

void inspiral_dist(double *deff, double *z, double *Vc, double m1_z,
		   double m2_z, double snr, double S_h[], int npoint,
		   double srate,double h100);

struct tm *utctime(const time_t *tp);
struct tm *gpstime(const time_t *tp);

/*
 *
 * MacOS support
 *
 */

/* Includes support for the MacOS  (c) S. Droz 1998
   You need libMac for this to work. */
#if __MACOS__ || macintosh
int GnuPlotCommand(char *command);
int PlayAudio(float *f, double Rate, int n); 
/* This is a kludge. the symbol macinosh will disapear and be replaced by
   __MACOS__. So let's get ready for the future. */
#ifndef __MACOS__ 
#define __MACOS__ 1
#endif
#endif 

#ifdef __cplusplus
}
#endif


/* close define of _GRASP_H at top of the header */
#endif

/*
 *
 *                  HEADER FOR THE TIMEFREQUENCY PACKAGE
 *
 *
 */

#ifndef PI
#define PI 3.1415926536
#endif

#define WIGNERTF 1       /* Wigner Transform with zero padding*/
#define WFFTWTF 2        /* Windowed Fourier Transform */
#define CHOIWILLIAMS 3   /* Choi William's transform */
#define WIGNERTF_NP 4       /* Wigner Transform With no zero padding*/

/* Structure definitions */
typedef struct {
    int run_number;        /*the directory number to store output files in  */
    float f_lower;         /* the lower frequency cutoff for the signals */
    int start_segment;     /* the segment to start analysing from */
    int transformtype;     /* the type of transform, WFT, WVD, CWD */
    int windowwidth;       /* the width of the window or sigma in case of the CWD */
    int offset_step_size;  /* used for averaging the time frequency map */
    int num_of_segments;   /* the number of segments to process */
    float maxpixelval;     /* computed from the value of RESCALE_FACTOR */
    int DIM;               /* the dimension of the segment */
    int ND;                /* the dimension of the subsegment */
    int PD;                /* the dimension of the time-freq map */
    int PRE;               /* to skip the initial PRESAFETY points in every segment */
    int POST;              /* to skip the final POSTSAFETY points in every segment */
    int TD;                /* = ND/PD */
    int FD;                /* = (ND/(4*PD)) */
    float rescale_factor;  /* rescale factor for rescaling images if desired */
    float hscale;          /* an arbitrary scaling number */
    int noisetype;         /*the type of noise */
    int srate;             /* the sampling rate in Hz.*/
}struct_tfparam;


typedef struct {
  double sigma;
  double low;
  double high;
  char   infile[256];
  char   outfile[256];
} dl_options;


/* Function prototypes */

void averageblock(float **, float *,struct_tfparam *);
void choiwill(unsigned long , float *, float *, int, struct_tfparam *);
float compute_scalefactor(float **,float, int);
void gen_quasiperiodic_signal(float  *, int ,float , float, float, float, float, float, int *);
int get_lines(double sigma, double high, double low, int rows, int cols, float **,float **);
int get_line_lens(double sigma, double high, double low, int rows, int cols, float **,char *);
void huetorgb(float,float *,float *,float *);
void normalize_picture(float **, int);
void pgmprint(float **,char *, int);
void plottf(float **,int);
void ppmprint(float **, char *, int);
void rescale(float **, int, float );
void time_freq_map(float *, struct_tfparam *,int,float **,float **);
void wigner(unsigned long , float *, float *, int);
void wigner_np(unsigned long , float *, float *, int);
void windowfft(long, float *, float *, int, struct_tfparam *);

    
/***************** END OF HEADER FOR TIMEFREQUENCY ROUTINES ********************/












