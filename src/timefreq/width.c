/*  Extract the width of a line for each line point; part of detect-lines.
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



#include "lines.h"



/* Constants */

/* Maximum search line length. */
#define MAX_LINE_WIDTH (2.5*sigma)

/* This constant is introduced because for very narrow lines the facet model
   width detection scheme sometimes extracts the line width too narrow.  Since
   the correction function has a very steep slope in that area, this will lead
   to lines of almost zero width, especially since the bilinear interpolation
   in correct.c will tend to overcorrect.  Therefore it is wise to make the
   extracted line width slightly larger before correction.  */
#define LINE_WIDTH_COMPENSATION 1.05

/* Minimum line width allowed (used for outlier check in fix_locations()) */
#define MIN_LINE_WIDTH 0.1

/* Maximum contrast allowed (used for outlier check in fix_locations()) */
#define MAX_CONTRAST 275.0



/* Local function prototypes */

static void fill_gaps __P((double *master, double *slave1, double *slave2,
                           contour *cont));
static void fix_locations __P((double *width_l, double *width_r,
                               double *grad_l, double *grad_r, double *pos_x,
                               double *pos_y, double *correction,
                               double *contr, double *asymm, double sigma,
                               long mode, bool correct_pos, contour *cont));



/* Modified Bresenham algorithm.  It returns in line all pixels that are
   intersected by a half line less than length away from the point (px,py)
   along the direction (nx,ny).  The point (px,py) must lie within the pixel
   of the origin, i.e., fabs(px) <= 0.5 and fabs(py) <= 0.5. */
void bresenham(nx, ny, px, py, length, line, num_points)
  double nx;
  double ny;
  double px;
  double py;
  double length;
  offset *line;
  long   *num_points;
{
  int i, n, x, y, s1, s2, xchg, maxit;
  double e, dx, dy, t;

  x = 0;
  y = 0;
  dx = ABS(nx);
  dy = ABS(ny);
  s1 = SGN(nx);
  s2 = SGN(ny);
  px *= s1;
  py *= s2;
  if (dy > dx) {
    t = dx;
    dx = dy;
    dy = t;
    t = px;
    px = py;
    py = t;
    xchg = 1;
  } else {
    xchg = 0;
  }
  maxit = ceil(length*dx);
  e = (0.5-px)*dy/dx-(0.5-py);
  n = 0;
  for (i=0; i<=maxit; i++) {
    line[n].x = x;
    line[n].y = y;
    n++;
    while (e >= -1e-8) {
      if (xchg) x += s1;
      else y += s2;
      e--;
      if (e > -1) {
        line[n].x = x;
        line[n].y = y;
        n++;
      }
    }
    if (xchg) y += s2;
    else x += s1;
    e += dy/dx;
  }
  *num_points = n;
}



/* Fill gaps in the arrays master, slave1, and slave2, i.e., points where
   master=0, by interpolation (interior points) or extrapolation (end points).
   The array master will usually be the width of the line, while slave1 and
   slave2 will be values that depend on master[i] being 0, e.g., the gradient
   at each line point.  The arrays slave1 and slave2 can be NULL. */
static void fill_gaps(master, slave1, slave2, cont)
  double  *master;
  double  *slave1;
  double  *slave2;
  contour *cont;
{
  long    i, j, k, s, e;
  long    num_points;
  double  m_s, m_e, s1_s, s1_e, s2_s, s2_e, d_r, d_c, arc_len, len;

  num_points = cont->num;
  for (i=0; i<num_points; i++) {
    if (master[i] == 0) {
      for (j=i+1; j<num_points; j++) {
        if (master[j] > 0)
          break;
      }
      m_s = 0;
      m_e = 0;
      s1_s = 0;
      s1_e = 0;
      s2_s = 0;
      s2_e = 0;
      if (i > 0 && j < num_points-1) {
        s = i;
        e = j-1;
        m_s = master[s-1];
        m_e = master[e+1];
        if (slave1 != NULL) {
          s1_s = slave1[s-1];
          s1_e = slave1[e+1];
        }
        if (slave2 != NULL) {
          s2_s = slave2[s-1];
          s2_e = slave2[e+1];
        }
      } else if (i > 0) {
        s = i;
        e = num_points-2;
        m_s = master[s-1];
        m_e = master[s-1];
        master[e+1] = m_e;
        if (slave1 != NULL) {
          s1_s = slave1[s-1];
          s1_e = slave1[s-1];
          slave1[e+1] = s1_e;
        }
        if (slave2 != NULL) {
          s2_s = slave2[s-1];
          s2_e = slave2[s-1];
          slave2[e+1] = s2_e;
        }
      } else if (j < num_points-1) {
        s = 1;
        e = j-1;
        m_s = master[e+1];
        m_e = master[e+1];
        master[s-1] = m_s;
        if (slave1 != NULL) {
          s1_s = slave1[e+1];
          s1_e = slave1[e+1];
          slave1[s-1] = s1_s;
        }
        if (slave2 != NULL) {
          s2_s = slave2[e+1];
          s2_e = slave2[e+1];
          slave2[s-1] = s2_s;
        }
      } else {
        s = 1;
        e = num_points-2;
        m_s = master[s-1];
        m_e = master[e+1];
        if (slave1 != NULL) {
          s1_s = slave1[s-1];
          s1_e = slave1[e+1];
        }
        if (slave2 != NULL) {
          s2_s = slave2[s-1];
          s2_e = slave2[e+1];
        }
      }
      arc_len = 0;
      for (k=s; k<=e+1; k++) {
        d_r = cont->row[k]-cont->row[k-1];
        d_c = cont->col[k]-cont->col[k-1];
        arc_len += sqrt(d_r*d_r+d_c*d_c);
      }
      len = 0;
      for (k=s; k<=e; k++) {
        d_r = cont->row[k]-cont->row[k-1];
        d_c = cont->col[k]-cont->col[k-1];
        len += sqrt(d_r*d_r+d_c*d_c);
        master[k] = (arc_len-len)/arc_len*m_s+len/arc_len*m_e;
        if (slave1 != NULL)
          slave1[k] = (arc_len-len)/arc_len*s1_s+len/arc_len*s1_e;
        if (slave2 != NULL)
          slave2[k] = (arc_len-len)/arc_len*s2_s+len/arc_len*s2_e;
      }
      i = j;
    }
  }
}



/* Correct the extracted line positions and widths.  The algorithm first closes
   gaps in the extracted data width_l, width_r, grad_l, and grad_r to provide
   meaningful input over the whole line.  Then the correction is calculated.
   After this, gaps that have been introduced by the width correction are again
   closed.  Finally, the position correction is applied if correct_pos is set.
   The results are returned in width_l, width_r, and cont. */
static void fix_locations(width_l,width_r,grad_l,grad_r,pos_x,pos_y,
                          correction,contr,asymm,sigma,mode,correct_pos,cont)
  double  *width_l;
  double  *width_r;
  double  *grad_l;
  double  *grad_r;
  double  *pos_x;
  double  *pos_y;
  double  *correction;
  double  *contr;
  double  *asymm;
  double  sigma;
  long    mode;
  bool    correct_pos;
  contour *cont;
{
  long    i;
  long    num_points;
  double  px, py;
  double  nx, ny;
  double  w_est, r_est, w_real, h_real, corr, w_strong, w_weak;
  double  correct, asymmetry, response, width, contrast;
  bool    weak_is_r;
  bool    correct_start, correct_end;

  fill_gaps(width_l,grad_l,NULL,cont);
  fill_gaps(width_r,grad_r,NULL,cont);

  num_points = cont->num;

  /* Calculate true line width, asymmetry, and position correction. */
  if (correct_pos) {
    /* Do not correct the position of a junction point if its width is found
       by interpolation, i.e., if the position could be corrected differently
       for each junction point, thereby destroying the junction. */
    correct_start = ((cont->cont_class == cont_no_junc ||
                      cont->cont_class == cont_end_junc ||
                      cont->cont_class == cont_closed) &&
                     (width_r[0] > 0 && width_l[0] > 0));
    correct_end = ((cont->cont_class == cont_no_junc ||
                    cont->cont_class == cont_start_junc ||
                    cont->cont_class == cont_closed) &&
                   (width_r[num_points-1] > 0 && width_l[num_points-1] > 0));
    /* Calculate the true width and assymetry, and its corresponding
       correction for each line point. */
    for (i=0; i<num_points; i++) {
      if (width_r[i] > 0 && width_l[i] > 0) {
        w_est = (width_r[i]+width_l[i])*LINE_WIDTH_COMPENSATION;
        if (grad_r[i] <= grad_l[i]) {
          r_est = grad_r[i]/grad_l[i];
          weak_is_r = TRUE;
        } else {
          r_est = grad_l[i]/grad_r[i];
          weak_is_r = FALSE;
        }
        line_corrections(sigma,w_est,r_est,&w_real,&h_real,&corr,
                         &w_strong,&w_weak);
        w_real /= LINE_WIDTH_COMPENSATION;
        corr /= LINE_WIDTH_COMPENSATION;
        width_r[i] = w_real;
        width_l[i] = w_real;
        if (weak_is_r) {
          asymm[i] = h_real;
          correction[i] = -corr;
        } else {
          asymm[i] = -h_real;
          correction[i] = corr;
        }
      }
    }

    fill_gaps(width_l,correction,asymm,cont);
    for (i=0; i<num_points; i++)
      width_r[i] = width_l[i];

    /* Adapt the correction for junction points if necessary. */
    if (!correct_start)
      correction[0] = 0;
    if (!correct_end)
      correction[num_points-1] = 0;

    for (i=0; i<num_points; i++) {
      px = pos_x[i];
      py = pos_y[i];
      nx = cos(cont->angle[i]);
      ny = sin(cont->angle[i]);
      px = px+correction[i]*nx;
      py = py+correction[i]*ny;
      pos_x[i] = px;
      pos_y[i] = py;
    }
  }

  /* Update the position of a line and add the extracted width. */
  cont->width_l = xcalloc(num_points,sizeof(float));
  cont->width_r = xcalloc(num_points,sizeof(float));
  for (i=0; i<num_points; i++) {
    cont->width_l[i] = width_l[i];
    cont->width_r[i] = width_r[i];
    cont->row[i] = pos_x[i];
    cont->col[i] = pos_y[i];
  }

  /* Now calculate the true contrast. */
  if (correct_pos) {
    cont->asymmetry = xcalloc(num_points,sizeof(float));
    cont->contrast = xcalloc(num_points,sizeof(float));
    for (i=0; i<num_points; i++) {
      response = cont->response[i];
      asymmetry = fabs(asymm[i]);
      correct = fabs(correction[i]);
      width = cont->width_l[i];
      if (width < MIN_LINE_WIDTH)
        contrast = 0;
      else
        contrast = 
          (response/fabs(phi2(correct+width,sigma)+
                         (asymmetry-1)*phi2(correct-width,sigma)));
      if (contrast > MAX_CONTRAST)
        contrast = 0;
      contr[i] = contrast;
    }
    fill_gaps(contr,NULL,NULL,cont);
    for (i=0; i<num_points; i++) {
      cont->asymmetry[i] = asymm[i];
      if (mode == MODE_LIGHT)
        cont->contrast[i] = contr[i];
      else
        cont->contrast[i] = -contr[i];
    }
  }
}



/* Extract the line width by using a facet model line detector on an image of
   the absolute value of the gradient. */
void compute_line_width(dx,dy,width,height,sigma,mode,correct_pos,contours,
                        num_contours)
  float   *dx;
  float   *dy;
  long    width;
  long    height;
  double  sigma;
  long    mode;
  bool    correct_pos;
  contour **contours;
  long    num_contours;
{
  float   *grad;
  long    i, j, k;
  long    r, c, l;
  long    x, y, dir;
  offset  *line;
  long    max_line, num_line;
  double  length;
  contour *cont;
  long    num_points, max_num_points;
  double  *width_r, *width_l;
  double  *grad_r, *grad_l;
  double  *pos_x, *pos_y, *correct, *asymm, *contrast;
  double  d, dr, dc, drr, drc, dcc;
  double  i1, i2, i3, i4, i5, i6, i7, i8, i9;
  double  t1, t2, t3, t4, t5, t6;
  double  eigval[2];
  double  eigvec[2][2];
  double  a, b, t;
  long    num;
  double  nx, ny;
  double  n1, n2;
  double  p1, p2;
  double  val;
  double  px, py;

  max_num_points = 0;
  for (i=0; i<num_contours; i++) {
    num_points = contours[i]->num;
    if (num_points > max_num_points)
      max_num_points = num_points;
  }

  width_l = xcalloc(max_num_points,sizeof(*width_l));
  width_r = xcalloc(max_num_points,sizeof(*width_r));
  grad_l = xcalloc(max_num_points,sizeof(*grad_l));
  grad_r = xcalloc(max_num_points,sizeof(*grad_r));
  pos_x = xcalloc(max_num_points,sizeof(*pos_x));
  pos_y = xcalloc(max_num_points,sizeof(*pos_y));
  correct = xcalloc(max_num_points,sizeof(*correct));
  contrast = xcalloc(max_num_points,sizeof(*contrast));
  asymm = xcalloc(max_num_points,sizeof(*asymm));

  grad = xcalloc(width*height,sizeof(*grad));
  memset(grad,0,width*height*sizeof(*grad));

  length = MAX_LINE_WIDTH;
  max_line = ceil(length*3);
  line = xcalloc(max_line,sizeof(*line));

  /* Compute the gradient image. */
  for (r=0; r<height; r++) {
    for (c=0; c<width; c++) {
      l = LINCOOR(r,c,width);
      grad[l] = sqrt(dx[l]*dx[l]+dy[l]*dy[l]);
    }
  }

  for (i=0; i<num_contours; i++) {
    cont = contours[i];
    num_points = cont->num;
    for (j=0; j<num_points; j++) {
      px = cont->row[j];
      py = cont->col[j];
      pos_x[j] = px;
      pos_y[j] = py;
      r = floor(px+0.5);
      c = floor(py+0.5);
      nx = cos(cont->angle[j]);
      ny = sin(cont->angle[j]);
      /* Compute the search line. */
      bresenham(nx,ny,0.0,0.0,length,line,&num_line);
      width_r[j] = width_l[j] = 0;
      /* Look on both sides of the line. */
      for (dir=-1; dir<=1; dir+=2) {
        for (k=0; k<num_line; k++) {
          x = BR(r+dir*line[k].x);
          y = BC(c+dir*line[k].y);
          i1 = grad[LINCOOR(BR(x-1),BC(y-1),width)];
          i2 = grad[LINCOOR(BR(x-1),y,width)];
          i3 = grad[LINCOOR(BR(x-1),BC(y+1),width)];
          i4 = grad[LINCOOR(x,BC(y-1),width)];
          i5 = grad[LINCOOR(x,y,width)];
          i6 = grad[LINCOOR(x,BC(y+1),width)];
          i7 = grad[LINCOOR(BR(x+1),BC(y-1),width)];
          i8 = grad[LINCOOR(BR(x+1),y,width)];
          i9 = grad[LINCOOR(BR(x+1),BC(y+1),width)];
          t1 = i1+i2+i3;
          t2 = i4+i5+i6;
          t3 = i7+i8+i9;
          t4 = i1+i4+i7;
          t5 = i2+i5+i8;
          t6 = i3+i6+i9;
          dr = (t3-t1)/6;
          dc = (t6-t4)/6;
          drr = (t1-2*t2+t3)/6;
          dcc = (t4-2*t5+t6)/6;
          drc = (i1-i3-i7+i9)/4;
          compute_eigenvals(2*drr,drc,2*dcc,eigval,eigvec);
          val = -eigval[0];
          if (val > 0.0) {
            n1 = eigvec[0][0];
            n2 = eigvec[0][1];
            a = 2.0*(drr*n1*n1+drc*n1*n2+dcc*n2*n2);
            b = dr*n1+dc*n2;
            solve_linear(a,b,&t,&num);
            if (num != 0) {
              p1 = t*n1;
              p2 = t*n2;
              if (fabs(p1) <= 0.5 && fabs(p2) <= 0.5) {
                /* Project the maximum point position perpendicularly onto the
                   search line. */
                a = 1;
                b = nx*(px-(r+dir*line[k].x+p1))+ny*(py-(c+dir*line[k].y+p2));
                solve_linear(a,b,&t,&num);
                d = (-i1+2*i2-i3+2*i4+5*i5+2*i6-i7+2*i8-i9)/9;
                if (dir == 1) {
                  grad_r[j] = d+p1*dr+p2*dc+p1*p1*drr+p1*p2*drc+p2*p2*dcc;
                  width_r[j] = fabs(t);
                } else {
                  grad_l[j] = d+p1*dr+p2*dc+p1*p1*drr+p1*p2*drc+p2*p2*dcc;
                  width_l[j] = fabs(t);
                }
                break;
              }
            }
          }
        }
      }
    }

    fix_locations(width_l,width_r,grad_l,grad_r,pos_x,pos_y,correct,contrast,
                  asymm,sigma,mode,correct_pos,cont);
  }

  free(line);
  free(grad);
  free(asymm);
  free(contrast);
  free(correct);
  free(pos_y);
  free(pos_x);
  free(grad_r);
  free(grad_l);
  free(width_r);
  free(width_l);
}
