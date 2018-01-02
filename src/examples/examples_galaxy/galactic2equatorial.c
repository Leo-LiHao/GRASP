/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"
#define DEG_TO_RAD  ((float)(M_PI/180))
#define RAD_TO_DEG  ((float)180/M_PI)
#define RAD_TO_HOUR ((float)12/M_PI)

int main(int argc, char *argv[])
{
  float l,b,alpha,delta;
  double hra,mra,sra,ddec,mdec,sdec;

  if (argc==3) {
    l = atof(argv[1]);
    b = atof(argv[2]);
  } else {
    fprintf(stderr,"usage: %s l b\n",argv[0]);
    fprintf(stderr,"  l:  Galactic longitude (degrees)\n");
    fprintf(stderr,"  b:  Galactic latitude (degrees)\n");
    return 1;
  }

  printf("l   (deg): %6.2f\n",l);
  printf("b   (deg): %+6.2f\n",b);

  l *= DEG_TO_RAD;
  b *= DEG_TO_RAD;

  galactic_to_equatorial(l,b,&alpha,&delta);

  alpha *= RAD_TO_HOUR;
  delta *= RAD_TO_DEG;

  sra = 60*modf(60*modf(alpha,&hra),&mra);
  sdec = 60*modf(fabs(60*modf(delta,&ddec)),&mdec);

  printf("RA  (hms):  %02u %02u %02u\n",(char)hra,(char)mra,(char)sra);
  printf("Dec (dms): %+02d %02u %02u\n",(int)ddec,(char)mdec,(char)sdec);

  return 0;
}
