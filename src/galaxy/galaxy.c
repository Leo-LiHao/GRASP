/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"

/* useful constants */
#define TWOPI 2*M_PI
#define PC_UNIT 1.0               /* the unit of one parsec (pc)             */
#define KPC_UNIT 1.0e3*PC_UNIT    /* the unit of one kiloparsec (kpc)        */
#define MPC_UNIT 1.0e6*PC_UNIT    /* the unit of one megaparsec (Mpc)        */
#define DEG_TO_RAD M_PI/180.0     /* conversion factor degrees -> radians    */
#define SEC_TO_RAD M_PI/43200.0   /* conversion factor seconds -> radians    */

/* Galactic coordinate constants */
#define R_SUN 8500.0*PC_UNIT    /* distance of sun from Galactic centre (pc) */
#define RA_NGP 192.25*DEG_TO_RAD  /* RA (1950) of North Galactic Pole (rad)  */
#define DEC_NGP 27.4*DEG_TO_RAD   /* Dec (1950) of North Galactic Pole (rad) */
#define L_ASCEND 33.0*DEG_TO_RAD  /* Galactic longitude of ascending node of
                                     Galactic plane on equator (1950) (rad)  */

/* NS-NS binary Galactic disk distribution parameters */
#define R_SCALE 4800.0*PC_UNIT    /* scale radius (pc)                       */
#define H_SCALE 1000.0*PC_UNIT    /* scale height (pc)                       */

/* computes the local sidereal time (decimal hours) at a given time
   at a given longitude (degrees West) */
float local_sidereal_time(time_t time, float longitude)
{
  long julday(int, int, int);

  struct tm *utc;
  double jd,ut,t,t0,gst,lst;

  /* convert longitude (degrees West) into decimal hours West */
  longitude *= 12.0/180.0;

  /* utc is coordinated universal (calendar) time
     ut is the universal time of day in decimal hours */
  utc = gmtime(&time);
  ut = (float)(utc->tm_sec)/3600.0;
  ut += (float)(utc->tm_min)/60.0;
  ut += (float)(utc->tm_hour);

  /* field tm_mon is months since Jan [0,11], julday requires month in [1,12];
     field tm_year is years since 1900, julday requires years since 1BC;
     we want the JD at 0h, which is the JD of noon of the previous calendar day
     plus one half (since JD starts at noon) */
  jd = (double)julday(utc->tm_mon + 1, utc->tm_mday - 1, utc->tm_year + 1900);
  jd += 0.5;

  /* compute Greenwitch Sidereal Time in decimal hours */
  t = (jd - 2451545.0)/36525.0;
  t0 = 6.697374558 + t*(2400.051336 + t*0.000025862);
  gst = t0 + 1.002737909*ut;

  /* compute Local Sidereal Time in decimal hours
     (negative sign since measuring longitude west) */
  lst = gst - longitude;

  /* put lst in range [0,24) */
  while (lst>=24.0) lst -= 24.0;
  while (lst<0.0) lst += 24.0;

  return (float)(lst);
}



/* produces a random (effective) amplitude and phase coefficients of a
   Galactic NS-NS binary at a given local sidereal time (seconds) for
   a detector at a given latitude (radians N) and with a given arm
   orientation (radians CCW from N) */
void mc_chirp(float time, float latitude, float orientation, long *seed,
	      float *invMpc, float *c0, float *c1)
{
  float ran1(long *);

  float u1,u2,u3,u4,u5;
  float rad,alpha,delta,alt,azi,theta,phi,psi,mu,plus,cross,x,y,z;

  /* generate the random numbers required */
  u1 = ran1(seed);
  u2 = ran1(seed);
  u3 = ran1(seed);
  u4 = ran1(seed);
  u5 = ran1(seed);

  /* compute the sky position of the binary and convert to horizon coords */
  (void)sky_position(u1,u2,u3,&rad,&alpha,&delta);
  equatorial_to_horizon(alpha,delta,time,latitude,&azi,&alt);

  /* convert horizon coords to detector response coords
     and compute the detector beam pattern response */
  theta = M_PI_2 - alt;
  phi = azi + orientation;
  psi = TWOPI*u4;            /* random polarization of source */
  beam_pattern(theta,phi,psi,&plus,&cross);

  /* compute the amplitude and phase of the chirp waveform */
  mu = 1 - 2*u5;             /* random cosine of inclination of binary */
  x = 0.5*plus*(1 + mu*mu);  /* factor of -2 omitted for optimal orientation */
  y = mu*cross;              /* factor of -2 omitted for optimal orientation */
  z = sqrt(x*x + y*y);
  *invMpc = MPC_UNIT*z/rad;  /* inverse distance (Mpc) of effective source */
  *c0 = -x/z;                /* phase coefficient 0 */
  *c1 = -y/z;                /* phase coefficient 1 */

  return;
}


/* computes the beam pattern functions, plus and cross, for some specified
   angles, theta, phi, and psi (in radians) */
void beam_pattern(float theta, float phi, float psi, float *plus, float *cross)
{
  float c=cos(theta),fac1=0.5*(1+c*c)*cos(2*phi),fac2=c*sin(2*phi);
  float c2psi=cos(2*psi),s2psi=sin(2*psi);

  *plus = fac1*c2psi - fac2*s2psi;
  *cross = fac1*s2psi + fac2*c2psi;
  return;
}


/* generates a random sky position of a NS-NS binary */
void sky_position(float u1, float u2, float u3,
		  float *rad, float *alpha, float *delta)
{
  float r,z,phi,rho2,l,b;

  r = R_SCALE*sqrt(-2*log(u1));
  /* need to worry about sign of z:
       if u2 > 0.5, use 2*(1-u2) and assume positve z;
       otherwise use 2*u2 and assume negative z. */
  if (u2>0.5) z = -H_SCALE*log(2*(1-u2));
  else z = H_SCALE*log(2*u2);
  phi = TWOPI*u3;

  rho2 = R_SUN*R_SUN + r*r - 2*R_SUN*r*cos(phi);
  *rad = sqrt(z*z + rho2);
  l = atan2(r*sin(phi),R_SUN - r*cos(phi));
  b = asin(z/(*rad));
  galactic_to_equatorial(l,b,alpha,delta);

  return;
}


/* converts from equatorial to horizon coordinates for a given
   time (sidereal seconds) and latitude (rad)
   note: alpha, delta, azi, and alt all in radians */
void equatorial_to_horizon(float alpha, float delta, float time, float lat,
			   float *azi, float *alt)
{
  float h = time*SEC_TO_RAD - alpha;  /* hour angle (rad) */

  *alt = asin(sin(delta)*sin(lat) + cos(delta)*cos(lat)*cos(h));
  *azi = atan2(-cos(delta)*cos(lat)*sin(h),sin(delta) - sin(lat)*sin(*alt));

  return;
}


/* converts from galactic coordinates to equatorial coordinates
   note: l, b, alpha, and delta all in radians */
void galactic_to_equatorial(float l, float b, float *alpha, float *delta)
{
  float lm=l-L_ASCEND;

  *delta = asin(cos(b)*cos(DEC_NGP)*sin(lm) + sin(b)*sin(DEC_NGP));
  *alpha = atan2(cos(b)*cos(lm),
		 sin(b)*cos(DEC_NGP) - cos(b)*sin(DEC_NGP)*sin(lm));
  *alpha += RA_NGP;

  return;
}
