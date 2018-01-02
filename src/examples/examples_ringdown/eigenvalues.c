/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"

int main(int argc,const char *argv[])
{
  float a,da=0.1,eigen[4];
  int s,l,m;

  /* process the command line arguments */
  if (argc==4) { /* correct number of arguments */
    s = atoi(argv[1]);
    l = atoi(argv[2]);
    m = atoi(argv[3]);
  } else { /* incorrect number of arguments */
    fprintf(stderr,"usage: qn_eigen_values s l m\n");
    return 1;
  }

  /* scan through the range of a */
  for (a=1-da;a>-1;a-=da) {
    /* compute the eigenvalues */
    if (a<0) {
      if (m==0) break;
      qn_eigenvalues(eigen,a,s,l,-m);
    } else {
      qn_eigenvalues(eigen,a,s,l,m);
    }
    /* print the eigenvalues */
    printf("%f\t%f\t%f\t%f\t%f\n",a,eigen[0],eigen[1],eigen[2],eigen[3]);
  }

  return 0;
}
