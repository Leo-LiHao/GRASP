/* GRASP: Copyright 1997,1998  Bruce Allen */
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>

static char *rcsid="$Id: test_transform.c,v 1.3 1998/02/25 23:44:24 ballen Exp $\n$Name: RELEASE_1_9_8 $";

int main(int argc, char **argv)
{
  float angle,c,s,co[7],coef[7];
  float a[2][2],g[2][2],g1[2][2],gamma[2][2][2],gamma1[2][2][2];
  int i,j,k,l,m,n;

  /* Get rotation angle and random number seed. */
  if(argc!=2){
    fprintf(stderr,"Usage: %s ANGLE\n%s\n",argv[0],rcsid);
    return 1;
  }
  angle=atof(argv[1]);
  c=cos(angle);
  s=sin(angle);

  /* Get initial coefficient vector. */
  co[0]=0.5435;
  co[1]=0.7538;
  co[2]=0.3469;
  co[3]=0.0985;
  co[4]=0.3746;
  co[5]=0.8346;
  co[6]=0.2344;

  fprintf(stdout,"Initial vector:\n");
  for(i=0;i<7;i++)
    fprintf(stdout,"%7.4f ",co[i]);
  fprintf(stdout,"\nTransformed vectors:\n");

  /* Transform by direct computation: */
  /* Transform quadratic components. */
  coef[0]=co[0]*c*c+co[1]*c*s+co[2]*s*s;
  coef[1]=2.0*(co[2]-co[0])*c*s+co[1]*(c*c-s*s);
  coef[2]=co[0]*s*s-co[1]*c*s+co[2]*c*c;

  /* Transform cubic components. */
  coef[3]=co[3]*c*c*c+co[4]*s*s*s+co[5]*c*c*s+co[6]*c*s*s;
  coef[4]=-co[3]*s*s*s+co[4]*c*c*c+co[5]*c*s*s-co[6]*c*c*s;
  coef[5]=(2.0*co[6]-3.0*co[3])*c*c*s+(3.0*co[4]-2.0*co[5])*c*s*s
    +co[5]*c*c*c-co[6]*s*s*s;
  coef[6]=(3.0*co[3]-2.0*co[6])*c*s*s+(3.0*co[4]-2.0*co[5])*c*c*s
    +co[5]*s*s*s+co[6]*c*c*c;

  for(i=0;i<7;i++)
    fprintf(stdout,"%7.4f ",coef[i]);
  fprintf(stdout,"\n");

  /* Transform by tensor formulation: */
  /* Compute tensors. */
  a[0][0]=a[1][1]=c;
  a[0][1]=-(a[1][0]=s);

  g[0][0]=co[0];
  g[0][1]=g[1][0]=co[1]/2.0;
  g[1][1]=co[2];

  gamma[0][0][0]=co[3];
  gamma[1][1][1]=co[4];
  gamma[0][0][1]=gamma[0][1][0]=gamma[1][0][0]=co[5]/3.0;
  gamma[0][1][1]=gamma[1][0][1]=gamma[1][1][0]=co[6]/3.0;

  /* Transform quadratic components. */
  for(i=0;i<2;i++)
    for(j=0;j<2;j++){
      g1[i][j]=0.0;
      for(l=0;l<2;l++)
	for(m=0;m<2;m++)
	  g1[i][j]+=g[l][m]*a[l][i]*a[m][j];
    }

  /* Transform cubic components. */
  for(i=0;i<2;i++)
    for(j=0;j<2;j++)
      for(k=0;k<2;k++){
	gamma1[i][j][k]=0.0;
	for(l=0;l<2;l++)
	  for(m=0;m<2;m++)
	    for(n=0;n<2;n++)
	      gamma1[i][j][k]+=gamma[l][m][n]*a[l][i]*a[m][j]*a[n][k];
      }

  /* Compute coefficient vector. */
  coef[0]=g1[0][0];
  coef[1]=g1[0][1]*2.0;
  coef[2]=g1[1][1];

  coef[3]=gamma1[0][0][0];
  coef[4]=gamma1[1][1][1];
  coef[5]=gamma1[0][0][1]*3.0;
  coef[6]=gamma1[0][1][1]*3.0;

  for(i=0;i<7;i++)
    fprintf(stdout,"%7.4f ",coef[i]);
  fprintf(stdout,"\n");

  return 0;
}
