/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"

static char *rcsid="$Id: detector_info.c,v 1.12 1999/09/08 22:55:03 ballen Exp $\
\n$Name: RELEASE_1_9_8 $";

void spline(float [],float [],int,float,float,float []);
void splint(float [],float [],float [],int,float,float *);

/* Function Name: detector_site() */

#define EQUATORIAL (6.37814e+08) /* earth's equatorial radius, in cm */
#define POLAR      (6.35676e+08) /* earth's polar radius, in cm */

void detector_site(char *detectors_file,int site_choice,
		   float site_parameters[9],char *site_name,
		   char *noise_file,char *whiten_file)
{
  int    i;
  int    site_id,neg;
  float  site_location_north,site_location_west,
         arm1_orientation,arm2_orientation,arm_length_cm;
  double factor,psi,phi,phi1,phi2,a2,b2,denom,
         location[3],east[3],north[3],arm1[3],arm2[3];

  char line[256];

  FILE *fp;
  fp=grasp_open("GRASP_PARAMETERS",detectors_file,"r");

  if ((neg=(site_choice<0))) site_choice = abs(site_choice);

  /* READ DETECTOR SITE DATA FROM FILE */

  while(1) {

    if (fgets(line,sizeof(line),fp)==NULL) {
	GR_start_error("detector_site()",rcsid,__FILE__,__LINE__);
	GR_report_error("%d is not a recognized site_id number.\n",site_choice);
	GR_end_error();
        abort();
    }
    if (line[0]!='#') {
      sscanf(line,"%d %e %e %e %e %e %s %s %s\n",
	     &site_id,&site_location_north,&site_location_west,
	     &arm1_orientation,&arm2_orientation,&arm_length_cm,
	     site_name,noise_file,whiten_file);
      if (site_choice==site_id) break;
    }
  }

  if (neg) { /* site_choice was negative: return raw parameters */
    *site_parameters++ = site_location_north;
    *site_parameters++ = site_location_west;
    *site_parameters++ = arm1_orientation;
    *site_parameters++ = arm2_orientation;
    *site_parameters   = arm_length_cm;
    return;
  }

  /* CALCULATE SITE PARAMETERS */

  factor = M_PI/180.0;
  psi=factor*site_location_north;
  phi= -1.0*factor*site_location_west;
  phi1=factor*arm1_orientation;
  phi2=factor*arm2_orientation;

  /* 
   Correct for oblateness of earth, use reference spheroid with   
   EQUATORIAL and POLAR radii given in table 4.2 of
   "Spacecraft attitude determination and control",  
   Ed. James R. Wertz, D. Reidel Publishing Co., Boston, 1978.
  */

  a2=EQUATORIAL*EQUATORIAL;
  b2=POLAR*POLAR;
  denom=sqrt(a2*cos(psi)*cos(psi)+b2*sin(psi)*sin(psi));

  location[0]=a2*cos(phi)*cos(psi)/denom;
  location[1]=a2*sin(phi)*cos(psi)/denom;
  location[2]=b2*sin(psi)/denom;

  /* unit length vectors in the north and east directions */
  north[0]= -1.0*cos(phi)*sin(psi);
  north[1]= -1.0*sin(phi)*sin(psi);
  north[2]= cos(psi);

  east[0]= -1.0*sin(phi);
  east[1]= cos(phi);
  east[2]= 0.0;

  for (i=0;i<3;i++) {
    arm1[i]=(arm_length_cm)*(cos(phi1)*north[i]-sin(phi1)*east[i]);
    arm2[i]=(arm_length_cm)*(cos(phi2)*north[i]-sin(phi2)*east[i]);
  }

  for (i=0;i<3;i++) {
    site_parameters[i]=location[i];
    site_parameters[i+3]=arm1[i];
    site_parameters[i+6]=arm2[i];
  }

  fclose(fp);

  return;
}
/* Function Name: noise_power() */

#define MAX 65536  /* (assumed) maximum size of data file */ 

void noise_power(char *noise_file,int n,float delta_f,double *power)
{
  int   i;
  int   file_length;
  float f,value,*x,*y,*y2;
  
  char line[256];

  FILE *fp;
  fp=grasp_open("GRASP_PARAMETERS",noise_file,"r");

  /* MEMORY ALLOCATION */
  x=(float *)malloc(MAX*sizeof(float));
  y=(float *)malloc(MAX*sizeof(float));
  y2=(float *)malloc(MAX*sizeof(float));
  if (! ( x && y && y2 ) ) {
    free(power);
    GR_start_error("noise_power()",rcsid,__FILE__,__LINE__);
    GR_report_error("Not enough memory available to read data.\n");
    GR_report_error("Freed memory pointed to by power!\n");
    GR_end_error();
    power=NULL;
    abort();
  }
  /* READ IN DATA FROM FILE AND STORE IN ARRAYS */
  i=0;

  while (1) {
    if (fgets(line,sizeof(line),fp)==NULL) {
      fclose(fp);
      break;
    }
    if (line[0]!='#') {
      sscanf(line,"%e %e\n",&x[i],&y[i]);
      i++;

      /* check that we have not overflowed arrays */
      if (i>=MAX) {
	GR_start_error("noise_power()",rcsid,__FILE__,__LINE__);
	GR_report_error("Internal array too small: MAX=%d.\n",MAX);
	GR_report_error("but file %s has more lines of data!\n",noise_file);
	GR_report_error("Recompile with larger value of MAX please.\n");
	GR_end_error();
        abort();
      }
    }
  }

  file_length=i;
	
  /* CALL SPLINE ROUTINE ONCE TO GET 2ND DERIVATIVES */
  spline(x-1,y-1,file_length,1.0e30,1.0e30,y2-1);

  
  /* CALL SPLINT ROUTINE FOR EACH DISCRETE FREQUENCY VALUE */
  for (i=0;i<n;i++) {

    f=i*delta_f;
    splint(x-1,y-1,y2-1,file_length,f,&value);

    /* square to get noise power in units of strain^2/Hz */
    /* (careful! about double vs float cast) */
    power[i]=value;
    power[i]*=power[i];
  }

  /* FREE STORAGE */
  free(x);
  free(y);
  free(y2);

  return;
}
/* Function Name: whiten() */

#define MAX 65536  /* (assumed) maximum size of data file */ 

void whiten(char *whiten_file,int n,float delta_f,double *whiten_out)
{
  int   i;
  int   file_length;
  float f,re_value,im_value;
  float *x,*re_y,*im_y,*re_y2,*im_y2;

  char line[256];

  FILE *fp;
  fp=grasp_open("GRASP_PARAMETERS",whiten_file,"r");

  /* MEMORY ALLOCATION */
  x=(float *)malloc(MAX*sizeof(float));
  re_y=(float *)malloc(MAX*sizeof(float));
  im_y=(float *)malloc(MAX*sizeof(float));
  re_y2=(float *)malloc(MAX*sizeof(float));
  im_y2=(float *)malloc(MAX*sizeof(float));

  /* READ IN DATA FROM FILE AND STORE IN ARRAYS */
  i=0;

  while (1) {
    if (fgets(line,sizeof(line),fp)==NULL) {
      fclose(fp);
      break;
    }
    if (line[0]!='#') {
      sscanf(line,"%e %e %e\n",&x[i],&re_y[i],&im_y[i]);
      i++;

      /* check that we have not overflowed arrays */
      if (i>=MAX) {
	GR_start_error("whiten()",rcsid,__FILE__,__LINE__);
	GR_report_error("Internal array too small: MAX=%d.\n",MAX);
	GR_report_error("but file %s has more lines of data!\n",whiten_file);
	GR_report_error("Recompile with larger value of MAX please.\n");
	GR_end_error();
        abort();
     }
    }
  }

  file_length=i;
	
  /* CALL SPLINE ROUTINE TO GET 2ND DERIVATIVES */
  spline(x-1,re_y-1,file_length,1.0e30,1.0e30,re_y2-1);
  spline(x-1,im_y-1,file_length,1.0e30,1.0e30,im_y2-1);

  
  /* CALL SPLINT ROUTINE FOR EACH DISCRETE FREQUENCY VALUE */
  for (i=0;i<n;i++) {

    f=i*delta_f;
    splint(x-1,re_y-1,re_y2-1,file_length,f,&re_value);
    splint(x-1,im_y-1,im_y2-1,file_length,f,&im_value);

    whiten_out[2*i]=re_value;
    whiten_out[2*i+1]=im_value;

  }

  /* FREE STORAGE */
  free(x);
  free(re_y);
  free(im_y);
  free(re_y2);
  free(im_y2);

  return;
}

/* Function Name: overlap() */
#define DOT(a,b) (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])

void overlap(float *site1_parameters,float *site2_parameters,int n,
	     float delta_f,double *gamma12)
{
  int    i;
  float  f;
  float  delta[3],x1[3],x2[3],y1[3],y2[3],m[3];
  float  x1x2,x1y2,y1y2,y1x2,x1m,x2m,y1m,y2m;
  float  c1,c2,c3,distance;
  double alpha,a2,a4,s,c;
  double bessel_j0,bessel_j1,bessel_j2;


  /* CALCULATE THE FREQUENCY INDEPENDENT COEFFICIENTS c1,c2,c3 */

  /* construct separation vector, arm vectors */
  for (i=0;i<3;i++) {
    delta[i]=site1_parameters[i]-site2_parameters[i];
    x1[i]=site1_parameters[3+i];
    y1[i]=site1_parameters[6+i];
    x2[i]=site2_parameters[3+i];
    y2[i]=site2_parameters[6+i];
  }
  
  /* construct UNIT VECTORS in the directions of the arms */
  normalize(x1);
  normalize(y1);
  normalize(x2);
  normalize(y2);

  /* unit separation vector between detectors */
  distance=sqrt(DOT(delta,delta));
  if (distance>0.0)
    for (i=0;i<3;i++) m[i]=delta[i]/distance;
  else 
    /* need separation vector in plane to define delta, Delta */
    for (i=0;i<3;i++) m[i]=x1[i];

  /* now construct the coefficients that will determine the */
  /* overlap reduction function */
  x1x2=DOT(x1,x2);
  y1y2=DOT(y1,y2);
  x1y2=DOT(x1,y2);
  y1x2=DOT(y1,x2);
  
  x1m=DOT(x1,m);
  x2m=DOT(x2,m);
  y1m=DOT(y1,m);
  y2m=DOT(y2,m);
  
  /* c1 = d1 : d2 */
  c1=x1x2*x1x2+y1y2*y1y2-x1y2*x1y2-y1x2*y1x2;
  c1*=0.25;

  /* c2 = m . d1 . d2 . m */
  c2=x1m*x2m*x1x2+y1m*y2m*y1y2-y1m*x2m*y1x2-x1m*y2m*x1y2;
  c2*=0.25;

  /* c3 = (m . d1 . m)(m . d2 . m) */
  c3= (x1m*x1m - y1m*y1m)*(x2m*x2m - y2m*y2m);
  c3*=0.25;

  if (distance==0.0) c2=c3=0.0;


  /* CALCULATE THE OVERLAP REDUCTION FUNCTION AT EACH DISCRETE FREQ */

  for (i=0;i<n;i++) {
    f=i*delta_f;
    alpha = 2*M_PI*f*distance/C_LIGHT;
    a2=alpha*alpha;

    /* if the argument is close to zero, use power series */
    if (alpha<0.01) {
      a4=a2*a2;
      bessel_j0=1.0-a2/6.0+a4/120.0;
      bessel_j1=(1.0-a2/10.0+a4/280.0)/3.0;
      bessel_j2=(1.0-a2/14.0+a4/504.0)/15.0;
    }
    else {
      s=sin(alpha);
      c=cos(alpha);

      /* the standard spherical bessel functions */
      bessel_j0=s/alpha;
      bessel_j1=(bessel_j0-c)/alpha;
      bessel_j2=(3.0*bessel_j1-s)/alpha;

      /* slightly modified */
      bessel_j1/=alpha;
      bessel_j2/=a2;
    }

    gamma12[i]=c1*(5.0*bessel_j0-10.0*bessel_j1+5.0*bessel_j2)+
               c2*(-10.0*bessel_j0+40.0*bessel_j1-50.0*bessel_j2)+
               c3*(2.5*bessel_j0-25.0*bessel_j1+87.5*bessel_j2);

  }
  return;
}	
/* Function Name: normalize() */

#define DOT(a,b) (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])

void normalize(float *v)
{
  float length;
  int i;
  
  length=DOT(v,v);
  length=1.0/sqrt(length);
  for (i=0;i<3;i++) v[i]*=length;

  return;
}
