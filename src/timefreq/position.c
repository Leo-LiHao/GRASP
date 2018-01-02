/*  Extract line points from an image; part of detect-lines.
    Copyright (C) 1996-1998 Carsten Steger

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2, or (at your option)
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA. */


/* Changes Made by R. Balasubramanian for incorporating the the detect lines code to incorporate
   within GRASP (May 10th 1999) 

   The type of the pointer image has been changed to float from unsigned char in all the function detect_lines()
  
   */



#include "lines.h"



/* Constants */

/* The pixel boundaries need to be enlarged slightly since in practice it
   frequently happens for neighboring pixels a and b that pixel a says a
   maximum lies within pixel b and vice versa.  This presents no problem since
   linking algoritm will take care of this. */
#define PIXEL_BOUNDARY 0.6



/* Local function prototypes */

static void compute_line_points __P((float *ku[5], byte *ismax, float *ev,
                                     float *nx, float *ny, float *px,
                                     float *py, long width, long height,
                                     double low, double high, long mode));



/* Solve the linear equation a*x+b=0 and return the result in t and the number
   of solutions in num. */
void solve_linear(a,b,t,num)
  double a;
  double b;
  double *t;
  long   *num;
{
  if (a == 0.0) {
    *num = 0;
    return;
  } else {
    *num = 1;
    *t = -b/a;
    return;
  }
}



/* Compute the eigenvalues and eigenvectors of the Hessian matrix given by
   dfdrr, dfdrc, and dfdcc, and sort them in descending order according to
   their absolute values. */
void compute_eigenvals(dfdrr,dfdrc,dfdcc,eigval,eigvec)
  double dfdrr;
  double dfdrc;
  double dfdcc;
  double eigval[2];
  double eigvec[2][2];
{
  double theta, t, c, s, e1, e2, n1, n2; /* , phi; */

  /* Compute the eigenvalues and eigenvectors of the Hessian matrix. */
  if (dfdrc != 0.0) {
    theta = 0.5*(dfdcc-dfdrr)/dfdrc;
    t = 1.0/(fabs(theta)+sqrt(theta*theta+1.0));
    if (theta < 0.0) t = -t;
    c = 1.0/sqrt(t*t+1.0);
    s = t*c;
    e1 = dfdrr-t*dfdrc;
    e2 = dfdcc+t*dfdrc;
  } else {
    c = 1.0;
    s = 0.0;
    e1 = dfdrr;
    e2 = dfdcc;
  }
  n1 = c;
  n2 = -s;

  /* If the absolute value of an eigenvalue is larger than the other, put that
     eigenvalue into first position.  If both are of equal absolute value, put
     the negative one first. */
  if (fabs(e1) > fabs(e2)) {
    eigval[0] = e1;
    eigval[1] = e2;
    eigvec[0][0] = n1;
    eigvec[0][1] = n2;
    eigvec[1][0] = -n2;
    eigvec[1][1] = n1;
  } else if (fabs(e1) < fabs(e2)) {
    eigval[0] = e2;
    eigval[1] = e1;
    eigvec[0][0] = -n2;
    eigvec[0][1] = n1;
    eigvec[1][0] = n1;
    eigvec[1][1] = n2;
  } else {
    if (e1 < e2) {
      eigval[0] = e1;
      eigval[1] = e2;
      eigvec[0][0] = n1;
      eigvec[0][1] = n2;
      eigvec[1][0] = -n2;
      eigvec[1][1] = n1;
    } else {
      eigval[0] = e2;
      eigval[1] = e1;
      eigvec[0][0] = -n2;
      eigvec[0][1] = n1;
      eigvec[1][0] = n1;
      eigvec[1][1] = n2;
    }
  }
}



/* For each point in the image determine whether there is a local maximum of
   the second directional derivative in the direction (nx[l],ny[l]) within the
   pixels's boundaries.  If so, set ismax[l] to 2 if the eigenvalue ev[l] is
   larger than high, to 1 if ev[l] is larger than low, and to 0 otherwise.
   Furthermore, put the sub-pixel position of the maximum into (px[l],py[l]).
   The parameter mode determines whether maxima (dark lines points) or minima
   (bright line points) should be selected.  The partial derivatives of the
   image are input as ku[]. */
static void compute_line_points(ku,ismax,ev,nx,ny,px,py,width,height,
                                low,high,mode)
  float   *ku[5];
  byte    *ismax;
  float   *ev;
  float   *nx;
  float   *ny;
  float   *px;
  float   *py;
  long    width;
  long    height;
  double  low;
  double  high;
  long    mode;
{
  long    r, c, l;
  double  k[5];
  double  eigval[2];
  double  eigvec[2][2];
  double  a, b, t;
  long    num;
  double  n1, n2;
  double  p1, p2;
  double  val;

  for (r=0; r<height; r++) {
    for (c=0; c<width; c++) {
      l = LINCOOR(r,c,width);
      k[0] = ku[0][l];
      k[1] = ku[1][l];
      k[2] = ku[2][l];
      k[3] = ku[3][l];
      k[4] = ku[4][l];
      ev[l] = 0.0;
      nx[l] = 0.0;
      ny[l] = 0.0;
      compute_eigenvals(k[2],k[3],k[4],eigval,eigvec);
      if (mode == MODE_LIGHT)
        val = -eigval[0];
      else
        val = eigval[0];
      if (val > 0.0) {
        ev[l] = val;
        n1 = eigvec[0][0];
        n2 = eigvec[0][1];
        a = k[2]*n1*n1+2.0*k[3]*n1*n2+k[4]*n2*n2;
        b = k[0]*n1+k[1]*n2;
        solve_linear(a,b,&t,&num);
        if (num != 0) {
          p1 = t*n1;
          p2 = t*n2;
          if (fabs(p1) <= PIXEL_BOUNDARY && fabs(p2) <= PIXEL_BOUNDARY) {
            if (val >= low) {
              if (val >= high)
                ismax[l] = 2;
              else
                ismax[l] = 1;
            }
            nx[l] = n1;
            ny[l] = n2;
            px[l] = r+p1;
            py[l] = c+p2;
          }
        }
      }
    }
  }
}



/* Main routine to detect lines in an image of dimension width * height.  The
   extracted lines are returned in result, while num_result is the number of
   detected lines.  The parameter sigma is the amount of smoothing that the
   Gaussian kernel performs, while low and high are the hysteresis thresholds
   used in the linking algorithm.  With mode, either bright or dark lines can
   be selected.  The parameter compute_width determines whether the line width
   should be extracted, while correct_pos determines whether the line width
   and position correction should be applied. */
void detect_lines(image,width,height,result,num_result,sigma,low,high,mode,
                  compute_width,correct_pos,extend_lines)
  float    *image;
  long    width;
  long    height;
  contour ***result;
  long    *num_result;
  double  sigma;
  double  low;
  double  high;
  long    mode;
  bool    compute_width;
  bool    correct_pos;
  bool    extend_lines;
{
  long   i;
  byte   *ismax;
  float  *ev, *n1, *n2, *p1, *p2;
  float  *k[5];

  for (i=0;i<5;i++)
    k[i] = xcalloc(width*height,sizeof(float));

  convolve_gauss(image,k[0],width,height,sigma,DERIV_R);
  convolve_gauss(image,k[1],width,height,sigma,DERIV_C);
  convolve_gauss(image,k[2],width,height,sigma,DERIV_RR);
  convolve_gauss(image,k[3],width,height,sigma,DERIV_RC);
  convolve_gauss(image,k[4],width,height,sigma,DERIV_CC);

  ismax = xcalloc(width*height,sizeof(*ismax));
  ev = xcalloc(width*height,sizeof(*ev));
  n1 = xcalloc(width*height,sizeof(*n1));
  n2 = xcalloc(width*height,sizeof(*n2));
  p1 = xcalloc(width*height,sizeof(*p1));
  p2 = xcalloc(width*height,sizeof(*p2));
  memset(ismax,0,width*height*sizeof(*ismax));
  memset(ev,0,width*height*sizeof(*ev));

  compute_line_points(k,ismax,ev,n1,n2,p1,p2,width,height,low,high,mode);
  compute_contours(ismax,ev,n1,n2,p1,p2,k[0],k[1],result,num_result,sigma,
                   extend_lines,mode,low,high,width,height);
  if (compute_width)
    compute_line_width(k[0],k[1],width,height,sigma,mode,correct_pos,*result,
                       *num_result);

  free(p2);
  free(p1);
  free(n2);
  free(n1);
  free(ev);
  free(ismax);

  for (i=4;i>=0;i--)
    free(k[i]);
}
