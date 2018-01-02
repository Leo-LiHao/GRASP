#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "grasp.h"

extern char *optarg;
extern int   optind;
int getopt(int, char * const [], const char *);

int main(int argc, char *argv[])
{
  void usage(char *);
  void realft(float [], unsigned long, int);
  void readNoisePower(FILE *, double [], int, float);
  void readCloseLimitData(FILE *, float [], int, float, float);
  void makeQuasiNormalRing(float [], int, float, float, float);

  const double hscale = 1e21;     /* a convenient scale factor            */
  const float  Msun   = 4.89e-6;  /* solar mass (s)                       */
  const float  srate  = 16384;    /* sample rate (Hz)                     */
  const int    npoint = 65536;    /* segment length (samples)             */

  double  scale;
  double *power;                  /* power spectrum %%$ S(f) $%%          */
  float  *weight;                 /* weigting factor %%$ 4/S(f) $%%       */
  float  *arrayA;                 /* waveform %%$ a(t) $%%                */
  float  *arrayB;                 /* waveform %%$ b(t) $%%                */
  float  *crossCorr;              /* %%$ \int df\,e^{2\pi ift}\tilde{a}(f)\tilde{b}^\ast(f)/S(f) $%% */ 
  float  *autoCorrA;              /* %%$ \int df\,e^{2\pi ift}\|\tilde{a}(f)\|^2/S(f) $%% */
  float  *autoCorrB;              /* %%$ \int df\,e^{2\pi ift}\|\tilde{b}(f)\|^2/S(f) $%% */
  float   spin = 0.35;            /* dimensionless spin of black hole     */
  float   mass = 200;             /* mass of black hole (solar masses)    */
  float   norm;
  float   max;
  char    powfile[256] = "ligo-0.dat";  /* noise power filename           */
  FILE   *fp;
  int     i;

  /* allocate memory to arrays */
  assert(power     = (double *)malloc((npoint/2+1)*sizeof(double)));
  assert(weight    = (float  *)malloc((npoint/2+1)*sizeof(float)));
  assert(arrayA    = (float  *)malloc(npoint*sizeof(float)));
  assert(arrayB    = (float  *)malloc(npoint*sizeof(float)));
  assert(crossCorr = (float  *)malloc(npoint*sizeof(float)));
  assert(autoCorrA = (float  *)malloc(npoint*sizeof(float)));
  assert(autoCorrB = (float  *)malloc(npoint*sizeof(float)));

  while (1) { /* parse command line options */

    int c;

    /* call the standard C library option parser */
    c = getopt(argc,argv,"hs:m:p:");
    if (c == -1) break;

    switch (c) {

    case 'h': /* print a simple message and exit */
      usage(argv[0]);
      exit(0);

    case 's':
      spin = atof(optarg);
      assert(spin >= 0 && spin < 1);
      break;

    case 'm':
      mass = atof(optarg);
      assert(mass > 0);
      break;

    case 'p':
      strncpy(powfile,optarg,sizeof(powfile));
      break;

    default:  /* something went wrong */
      fprintf(stderr,"warning: getopt returned character code O%o\n",c);

    }

  }

  switch (argc - optind) { /* process remaining command line arguments */

  case 0: /* compare a default data file and a computed ringdown */
    fprintf(stderr,"compare waveform data from file close-limit.dat\n");
    assert(fp = fopen("close-limit-insp.dat","r"));
    readCloseLimitData(fp,arrayA,npoint,srate,mass*Msun);
    fclose(fp);
    fprintf(stderr,"with a computed ringdown waveform\n");
    makeQuasiNormalRing(arrayB,npoint,srate,mass*Msun,spin);
    fprintf(stderr,"  black hole mass: %.2f solar masses\n",mass);
    fprintf(stderr,"  black hole spin: %05.2f%% of extreme\n",spin*100);
    break;

  case 1: /* compare a specified data file and a computed ringdown */
    fprintf(stderr,"compare waveform data from file %s\n",argv[optind]);
    assert(fp = fopen(argv[optind++],"r"));
    readCloseLimitData(fp,arrayA,npoint,srate,mass*Msun);
    fclose(fp);
    fprintf(stderr,"with a computed ringdown waveform\n");
    makeQuasiNormalRing(arrayB,npoint,srate,mass*Msun,spin);
    fprintf(stderr,"  black hole mass: %.2f solar masses\n",mass);
    fprintf(stderr,"  black hole spin: %05.2f%% of extreme\n",spin*100);
    break;

  case 2: /* compare two specified data files */
    fprintf(stderr,"compare waveform data from file %s\n",argv[optind]);
    assert(fp = fopen(argv[optind++],"r"));
    readCloseLimitData(fp,arrayA,npoint,srate,mass*Msun);
    fclose(fp);
    fprintf(stderr,"with waveform data from file %s\n",argv[optind]);
    assert(fp = fopen(argv[optind++],"r"));
    readCloseLimitData(fp,arrayB,npoint,srate,mass*Msun);
    fprintf(stderr,"  black hole mass: %.2f solar masses\n",mass);
    fclose(fp);
    break;

  default: /* too many arguments */
    usage(argv[0]);
    exit(1);

  }

  /* get correlation weighting factor */
  fprintf(stderr,"using power spectrum from file %s\n",powfile);
  assert(fp = grasp_open("GRASP_PARAMETERS",powfile,"r"));
  readNoisePower(fp,power,npoint/2,srate);
  fclose(fp);
  scale = 4/(npoint*srate*hscale*hscale);
  weight[0] = weight[npoint/2] = 0;
  for (i = 1; i < npoint/2; ++i)
    weight[i] = scale/power[i];

  /* FFT waveform arrays */
  realft(arrayA-1,npoint,1);
  realft(arrayB-1,npoint,1);

  /* compute cross- and auto-correlations */
  autoCorrA[0] = autoCorrA[1] = 0;
  autoCorrB[0] = autoCorrB[1] = 0;
  crossCorr[0] = crossCorr[1] = 0;
  for (i = 1; i < npoint/2; ++i) {
    int ir = i + i;
    int ii = ir + 1;
    float ar = arrayA[ir];
    float ai = arrayA[ii];
    float br = arrayB[ir];
    float bi = arrayB[ii];
    float fac = weight[i];
    autoCorrA[ir] = fac*(ar*ar + ai*ai);
    autoCorrA[ii] = 0;
    autoCorrB[ir] = fac*(br*br + bi*bi);
    autoCorrB[ii] = 0;
    crossCorr[ir] = fac*(ar*br + ai*bi);
    crossCorr[ii] = fac*(ai*br - ar*bi);
  }
  realft(autoCorrA-1,npoint,-1);
  realft(autoCorrB-1,npoint,-1);
  realft(crossCorr-1,npoint,-1);

  /* compute fitting factor normalization */
  assert(autoCorrA[0] > 0);
  assert(autoCorrB[0] > 0);
  norm = sqrt(autoCorrA[0]*autoCorrB[0]);
  assert(norm > 0);

  /* find maximum of cross-correlation and print fitting factor */
  max = 0;
  for (i = 0; i < npoint; ++i)
    if (fabs(crossCorr[i]) > fabs(max))
      max = crossCorr[i];
  printf("fitting factor: %.2f%%\n",100*max/norm);

  return 0;
}


/* Print a message describing the usage of this program.
 * Arguments:
 *   @ *program @ the program name
 */
void usage(char *program)
{
  fprintf(stderr,"usage: %s [options] [file1 [file2]]\n",program);
  fprintf(stderr,"options:\n");
  fprintf(stderr,"  -h       prints this message\n");
  fprintf(stderr,"  -s spin  dimensionless spin of the black hole [0,1)\n");
  fprintf(stderr,"  -m mass  mass of the black hole (solar masses)\n");
  fprintf(stderr,"  -p file  file containing the noise power spectrum\n");
  return;
}


/* Read a data file containing a close-limit waveform.
 * Arguments:
 *   @ *fp     @ the data file
 *   @ *arr    @ array to store the data
 *   @  npoint @ size of array @*arr@
 *   @  srate  @ sample rate for data in @*arr@ (in Hz)
 *   @  mass   @ mass of black hole (in seconds)
 */
void readCloseLimitData(FILE *fp, float *arr, int npoint, float srate,
                       float mass)
{
  void spline(float [], float [], int, float, float, float []);
  void splint(float [], float [], float [], int, float, float *);

  const int   ninc    = 1024;
  const float natural = 1e30;
  float *time, *data, *datapp;
  int i, imax, nmax = ninc, n = 0;

  /* alocate memory */
  assert(data = (float *)malloc(nmax*sizeof(float)));
  assert(time = (float *)malloc(nmax*sizeof(float)));

  /* read waveform data file */
  while (EOF != fscanf(fp,"%e\t%e\n",time+n,data+n))
    if (++n >= nmax) {
      nmax += ninc;
      assert(data = (float *)realloc(data,nmax*sizeof(float)));
      assert(time = (float *)realloc(time,nmax*sizeof(float)));
    }

  /* use cubic spline interpolation to generate waveform at
   * the required sample times
   */
  assert(datapp = (float *)malloc(n*sizeof(float)));
  spline(time-1,data-1,n,natural,natural,datapp-1);
  imax = (int)floor((time[n-1] - time[0])*srate*mass);
  for (i = 0; i < npoint; ++i)
    if (i > imax)
      arr[i] = 0;
    else
      splint(time-1,data-1,datapp-1,n,(float)i/(srate*mass)+time[0],arr+i);

  /* free memory and return */
  free(time);
  free(data);
  free(datapp);
  return;
}

/* Generate a quasinormal ringdown waveform.
 * Arguments:
 *   @ *arr    @ array to store the data
 *   @  npoint @ size of array @*arr@
 *   @  srate  @ sample rate for data in @*arr@ (in Hz)
 *   @  mass   @ mass of black hole (in seconds)
 *   @  spin   @ dimensionless spin of black hole
 */
void makeQuasiNormalRing(float *arr, int npoint, float srate,
                         float mass, float spin)
{
  const float pi   = 3.14159265358979;
  const float freq = (1 - 0.63*pow(1-spin,0.3))/(2*pi*mass);
  const float qual = 2*pow(1-spin,-0.45);
  /* @freq@ and @qual@ are computed using the fits by Echeverria */
  int i;

  for (i = 0; i < npoint; ++i) {
    float time = (float)i/srate;
    arr[i] = exp(-pi*time*freq/qual)*cos(2*pi*time*freq);
  }

  return;
}

/* Read the noise power spectrum.
 * Arguments:
 *   @ *fp    @ the noise power file
 *   @ *power @ array to store the data
 *   @  n     @ size of array @*power@
 *   @  srate @ sample rate for data (in Hz)
 */
void readNoisePower(FILE *fp, double *power, int n, float srate)
{
  void spline(float [], float [], int, float, float, float []);
  void splint(float [], float [], float [], int, float, float *);
  const int   nmax    = 65536; /* assumed maximum size of data file */
  const float natural = 1e30;
  float *freq, *ampl, *amplpp;
  char line[100];
  int i, length;

  /* allocate memory */
  assert(freq   = (float *)malloc(nmax*sizeof(float)));
  assert(ampl   = (float *)malloc(nmax*sizeof(float)));
  assert(amplpp = (float *)malloc(nmax*sizeof(float)));

  /* read power spectrum data file */
  length = 0;
  while (1) {
    if (fgets(line,sizeof(line),fp) == NULL) /* end of file */
      break;
    if (line[0] != '#') {
      assert(length < nmax);
      sscanf(line,"%e\t%e\n",freq+length,ampl+length);
      ++length;
    }
  }

  /* use cubic spline interpolation to get the spectrum
   * at the required frequencies
   */
  spline(freq-1,ampl-1,length,natural,natural,amplpp-1);
  for (i = 0; i < n; ++i) {
    float f = i*srate/(float)n;
    float value;
    double dvalue;
    splint(freq-1,ampl-1,amplpp-1,length,f,&value);
    dvalue = (double)value;
    power[i] = dvalue*dvalue;
  }

  /* free memory and return */
  free(freq);
  free(ampl);
  free(amplpp);
  return;
}

/* The following routines is to be run to display data during
   debugging it's based on the GRASP graph() routine */

void graphSpec(float arr[], int n)
{
  FILE *fp;
  int i;
  fp = fopen("temp.graph","w");
  for (i = 0; i < n/2; ++i) {
    int ir = i + i;
    int ii = ir + 1;
    float re = arr[ir];
    float im = arr[ii];
    float pow = re*re + im*im;
    float arg = atan2(im,re);
    fprintf(fp,"%d\t%e\t%e\n",i,pow,arg);
  }
  fclose(fp);
  system("xmgr -nxy temp.graph 1>/dev/null 2>&1 &");
  return;
}
