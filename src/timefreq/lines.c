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

   The type of the pointer image has been changed to float from unsigned char in function 
   get_lines();

   new function added by Warren Anderson  get_line_lens()
   */


#include "lines.h"



/* Type declarations */

/* Command line options for the program */
typedef struct {
  double sigma;
  double low;
  double high;
  long   mode;
  bool   correct;
  bool   width;
  bool   extend;
  bool   postscript;
  bool   encapsulated;
  bool   image;
} options;

/* Local function prototypes */

static void check_sigma __P((double sigma, long width, long height));

/* Check whether sigma is in the allowed range for the current image size. */
static void check_sigma(sigma, width, height)
  double sigma;
  long width;
  long height;
{
  long min_dim;

  min_dim = width < height ? width : height;
  if (sigma < 0.4)
    print_error(ERR_SOR,"< 0.4");
  if (MASK_SIZE(MAX_SIZE_MASK_2,sigma) >= min_dim)
    print_error(ERR_SOR,"too large for image size");
}



/* The main program */
int get_lines(double sigma, double high, double low, int rows,
		int cols, float  **in_img, float **out_img){
  
  float    *image;
  contour **contours, *cont;
  long    i, num_cont;
  static options opts = {
    -1.0, -1.0, -1.0, MODE_LIGHT, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE
  };
  int retvalue;
  

  opts.sigma=sigma;
  opts.high=high;
  opts.low=low;

  check_sigma(opts.sigma,cols,rows);
  {
      int i,j;
      image = (float *) malloc(rows*cols*sizeof(float));
      for (i=0;i<rows;i++)
          for(j=0;j<cols;j++)
              *(image + i*cols + j) = *(in_img[i] + j);
  }
  
  detect_lines(image,cols,rows,&contours,&num_cont,opts.sigma,opts.low,
               opts.high,opts.mode,opts.width,opts.correct,opts.extend);
  retvalue = num_cont;
  
  show_lines(rows,cols,contours,num_cont,out_img);

  for (i=0; i<num_cont; i++) {
    cont = contours[i];
    if (cont->width_r != NULL)
      free(cont->width_r);
    if (cont->width_l != NULL)
      free(cont->width_l);
    if (cont->angle != NULL)
      free(cont->angle);
    if (cont->col != NULL)
      free(cont->col);
    if (cont->row != NULL)
      free(cont->row);
    free(cont);
  }
  free(contours);
  free(image);
  return retvalue;
}


int get_line_lens(double sigma, double high, double low, int rows,
		int cols, float **in_img, char *outfilename){
    float *image;
  contour **contours, *cont;
  long    i, num_cont;
  static options opts = {
    -1.0, -1.0, -1.0, MODE_LIGHT, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE
  };
  int retvalue;


  opts.sigma=sigma;
  opts.high=high;
  opts.low=low;

  check_sigma(opts.sigma,cols,rows);
  {
      int i,j;
      image = (float *) malloc(rows*cols*sizeof(float));
      for (i=0;i<rows;i++)
          for(j=0;j<cols;j++)
              *(image + i*cols + j) =  *(in_img[i] + j);
  }
    

  detect_lines(image,cols,rows,&contours,&num_cont,opts.sigma,opts.low,
               opts.high,opts.mode,opts.width,opts.correct,opts.extend);
  {
      FILE *fout;
      int j,num_points,intx,inty;
      float strength;
      
      if(num_cont!= 0){
          fout=fopen(outfilename,"w");
          fprintf(fout,"%ld\n",num_cont);
          for(i=0;i<num_cont;i++){
              cont = contours[i];
              num_points = cont->num;
              strength = 0.0;
              for(j=0;j<num_points;j++){
                  intx = (int) (cont->row[j] +0.5);
                  inty = (int) (cont->col[j] +0.5);
                  strength += *(in_img[intx] + inty);
              }
              fprintf(fout,"%d %f \n",num_points,strength);
          }
          fclose(fout);
      }
      
  }
  retvalue = num_cont;
  
  for (i=0; i<num_cont; i++) {
    cont = contours[i];
    if (cont->width_r != NULL)
      free(cont->width_r);
    if (cont->width_l != NULL)
      free(cont->width_l);
    if (cont->angle != NULL)
      free(cont->angle);
    if (cont->col != NULL)
      free(cont->col);
    if (cont->row != NULL)
      free(cont->row);
    free(cont);
  }
  free(contours);
  free(image);
  return retvalue;
}

extern void show_lines(int rows, int cols, contour **contours, long num_cont, 
                       float **out_img)
{
    int i,j,points;
    contour *cont;
    
    
    for(i=0;i<num_cont;i++){
		cont=contours[i];
        points=cont->num;
        for(j=0;j<points;j++){
            out_img[(int)(cont->row[j]+0.5)][(int)(cont->col[j]+0.5)]=1.0;
        }
    }
}


