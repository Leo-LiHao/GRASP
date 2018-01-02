/*  Threshold an image; part of detect-lines.
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



/* Threshold an image above min and return the result as a run-length encoded
   region in out. */
void threshold(image,min,width,height,out)
  byte   *image;
  long   min;
  long   width;
  long   height;
  region *out;
{
  long   grey;
  long   r,c,l,num,num_max;
  bool   inside;
  chord  *rl;

  inside = FALSE;
  num = 0;
  num_max = INITIAL_SIZE;
  rl = xcalloc(num_max,sizeof(chord));

  out->rl = NULL;
  out->num = 0;

  for (r=0; r<height; r++) {
    for (c=0; c<width; c++) {
      l = LINCOOR(r,c,width);
      grey = image[l];
      if (grey >= min) {
        if (!inside) {
          inside = TRUE;
          rl[num].r = r;
          rl[num].cb = c;
        }
      } else {
        if (inside) {
          inside = FALSE;
          rl[num].ce = c - 1;
          num++;
          if (num >= num_max) {
            num_max = floor((double)(num_max*REALLOC_FACTOR));
            rl = xrealloc(rl,num_max*sizeof(chord));
          }
        }
      }
    }
    if (inside) {
      inside = FALSE;
      rl[num].ce = width-1;
      num++;
      if (num >= num_max) {
        num_max = floor((double)(num_max*REALLOC_FACTOR));
        rl = xrealloc(rl,num_max*sizeof(chord));
      }
    }
  }
  out->rl = xcalloc(num,sizeof(chord));
  memcpy(out->rl,rl,num*sizeof(chord));
  out->num = num;

  free(rl);
}
