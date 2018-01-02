/*  detect-lines, extract lines and their width from images.
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

   The type of the pointer image has been changed to float from unsigned char in all the relevant
   function prototypes.


   */


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#ifdef STDC_HEADERS
#include <stdlib.h>
#include <string.h>
#include <float.h>
#else /* not STDC_HEADERS */
#ifdef HAVE_STRING_H
#include <string.h>
#else /* not HAVE_STRING_H */
#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif /* HAVE_STRINGS_H */
#endif /* not HAVE_STRING_H */
#ifdef HAVE_MEMORY_H
#include <memory.h>
#endif /* HAVE_MEMORY_H */
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif /* HAVE_MALLOC_H */
#ifndef HAVE_STRCHR
#define strchr  index
#define strrchr rindex
#endif /* not HAVE_STRCHR */
extern char *strchr(), *strrchr();
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif /* not HAVE_MEMCPY */
#ifndef HAVE_MEMSET
#define memset(s,c,n) bzero((s),(n))
#endif /* not HAVE_MEMSET */
#endif /* not STDC_HEADERS */
#ifndef DBL_MAX
#define DBL_MAX 1.7976931348623157E+308
#endif /* not DBL_MAX */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif /* not M_PI */
#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880
#endif /* not M_SQRT2 */

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#include <pwd.h>
#endif

#include <errno.h>
/*extern char *sys_errlist[];*/

#ifndef PI
#define PI M_PI
#endif /* not PI */
#ifndef SQRT2
#define SQRT2 M_SQRT2
#endif /* not SQRT2 */

#ifdef __STDC__
#ifdef __P
#undef __P
#endif
#define __P(args) args
#else /* not __STDC__ */
#define __P(args) ()
#endif /* not __STDC__ */



/* Constants */

#define VERSION  1  /* Program version */
#define REVISION 2  /* Program revision */

#define FALSE 0
#define TRUE  1

#define DERIV_R  1  /* Derivative in row direction */
#define DERIV_C  2  /* Derivative in column direction */
#define DERIV_RR 3  /* Second derivative in row direction */
#define DERIV_RC 4  /* Second derivative in row and column direction */
#define DERIV_CC 5  /* Second derivative in column direction */

#define MODE_LIGHT 1  /* Extract bright lines */
#define MODE_DARK  2  /* Extract dark lines */

#define INITIAL_SIZE   100  /* Initial size of dynamically realloc'ed arrays */
#define REALLOC_FACTOR 2    /* Size increment for dynamic arrays */


/* Mask sizes in convol.c */

#define MAX_SIZE_MASK_0  3.09023230616781    /* Size for Gaussian mask */
#define MAX_SIZE_MASK_1  3.46087178201605    /* Size for 1st derivative mask */
#define MAX_SIZE_MASK_2  3.82922419517181    /* Size for 2nd derivative mask */

#define MASK_SIZE(MAX,sigma) ceil(MAX*sigma) /* Maximum mask index */


/* Error messages */

#define ERR_NOMEM "Out of memory"
#define ERR_FNF   "Could not open file"
#define ERR_NOPGM "Not a valid pgm file:"
#define ERR_SNS   "Sigma not specified"
#define ERR_SOR   "Sigma out of range:"
#define ERR_LNS   "Low not specified"
#define ERR_LOR   "Low out of range:"
#define ERR_HNS   "High not specified"
#define ERR_HOR   "High out of range:"
#define ERR_LGH   "Low > High"
#define ERR_CNW   "Line position correction impossible without line width"
#define ERR_INP   "Include-image option requires PostScript output"
#define ERR_UKO   "Unknown option:"
#define ERR_TMF   "Too many files specified:"



/* Macros */

/* Translate row and column coordinates of an image into an index into its
   one-dimensional array. */
#define LINCOOR(row,col,width) (long)((row)*(width)+(col))

/* Mirror the row coordinate at the borders of the image; height must be a
   defined variable in the calling function containing the image height. */
#define BR(row) ((row) < 0 ? -(row) : \
                 (row) >= height ? height - (row) + height - 2 : (row))

/* Mirror the column coordinate at the borders of the image; width must be a
   defined variable in the calling function containing the image width. */
#define BC(col) ((col) < 0 ? -(col) : \
                 (col) >= width ? width - (col) + width - 2 : (col))

/* Absolute value of x */
#ifndef ABS
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#endif

/* Sign of x */
#ifndef SGN
#define SGN(x) ((x) == 0 ? 0 : ((x) > 0 ? 1 : -1))
#endif



/* Type declarations */

/* Type for images: images are declared as byte* */
typedef unsigned char byte;

/* Boolean type */
typedef int bool;

/* Offsets to a specific location in the image.  An array of this type is
   returned by the modified Bresenham algorithm in width.c.  It is also used
   in link.c to hold an array of pixel locations to check for appropriate
   neighbors. */
typedef struct {
  long x;
  long y;
} offset;

/* This type determines the class of a line, i.e., whether its end points are
   junction points or whether the line forms a closed loop. */
typedef enum {
  cont_no_junc,    /* no end point is a junction */
  cont_start_junc, /* only the start point of the line is a junction */
  cont_end_junc,   /* only the end point of the line is a junction */
  cont_both_junc,  /* both end points of the line are junctions */
  cont_closed      /* the contour is closed */
} contour_class;

/* This type holds one extracted line.  The field num contains the number of
   points in the line.  The coordinates of the line points are given in the
   arrays row and col.  The array angle contains the direction of the normal
   to each line point, as measured from the row-axis.  Some people like to
   call the col-axis the x-axis and the row-axis the y-axis, and measure the
   angle from the x-axis.  To convert the angle into this convention, subtract
   PI/2 from the angle and normalize it to be in the interval [0,2*PI).  The
   array response contains the response of the operator, i.e., the second
   directional derivative in the direction of angle, at each line point.  The
   arrays width_l and width_r contain the width information for each line point
   if the algorithm was requested to extract it; otherwise they are NULL.  If
   the line position and width correction was applied the contents of width_l
   and width_r will be identical.  The arrays asymmetry and contrast contain
   the true asymmetry and contrast of each line point if the algorithm was
   instructed to apply the width and position correction.  Otherwise, they are
   set to NULL.  If the asymmetry, i.e., the weaker gradient, is on the right
   side of the line, the asymmetry is set to a positive value, while if it is
   on the left side it is set to a negative value. */
typedef struct {
  long  num;                /* number of points */
  float *row;               /* row coordinates of the line points */
  float *col;               /* column coordinates of the line points */
  float *angle;             /* angle of normal (measured from the row axis) */
  float *response;          /* response of line point (second derivative) */
  float *width_l;           /* width to the left of the line */
  float *width_r;           /* width to the right of the line */
  float *asymmetry;         /* asymmetry of the line point */
  float *contrast;          /* contrast of the line point */
  contour_class cont_class; /* contour class (e.g., closed, no_junc) */
} contour;

/* A chord in a run-length encoded region */
typedef struct { 
  short r;   /* row coordinate of the chord */
  short cb;  /* column coordinate of the start of the chord */
  short ce;  /* column coordinate of the end of the chord */
} chord;

/* Run-length encoded region of an image.  This type is returned by the
   threshold() function.  It provides the means to efficiently link line points
   into lines. */
typedef struct { 
  long  num;      /* number of chords */     
  chord *rl;      /* array of chords */
} region;



/* Global variables */

/* The basename of the program as called from the shell */
extern char *program_name;



/* Global function prototypes */

void *xmalloc __P((size_t size));
void *xcalloc __P((size_t nelem, size_t elsize));
void *xrealloc __P((void *ptr, size_t size));
extern double normal __P((double x));
extern void solve_linear __P((double a, double b, double *t, long *num));
extern void compute_eigenvals __P((double dfdrr, double dfdrc, double dfdcc,
                                   double eigval[2], double eigvec[2][2]));
extern void threshold __P((byte *image, long min, long width, long height,
                           region *out));
extern void bresenham __P((double nx, double ny, double px, double py,
                           double length, offset *line, long *num_points));
extern double phi0 __P((double x, double sigma));
extern double phi1 __P((double x, double sigma));
extern double phi2 __P((double x, double sigma));
extern void convolve_gauss __P((float *image, float *k, long width, long height,
                                double sigma, long deriv_type));
extern bool line_corrections __P((double sigma, double w_est, double r_est,
                                  double *w, double *h,  double *correction,
                                  double *w_strong, double *w_weak));
extern void compute_contours __P((byte *ismax, float *eigval, float *normx,
                                  float *normy, float *posx, float *posy,
                                  float *gradx, float *grady,
                                  contour ***result, long *num_result,
                                  double sigma, bool extend_lines, int mode,
                                  double low, double high, long width,
                                  long height));
extern void compute_line_width __P((float *dx, float *dy, long width,
                                    long height, double sigma, long mode,
                                    bool correct_pos, contour **contours,
                                    long num_contours));
extern void detect_lines __P((float *image, long width, long height,
                              contour ***result, long *num_result,
                              double sigma, double low, double high, long mode,
                              bool compute_width, bool correct_pos,
                              bool extend_lines));
extern int conv_img(int rows, int cols, byte **in_img, byte *image);
extern void show_lines(int rows, int cols, contour **contours, long num_cont, 
				float **out_img);
extern void print_error __P((const char *message, const char *info));
extern void print_ascii( contour **, long , char *);






