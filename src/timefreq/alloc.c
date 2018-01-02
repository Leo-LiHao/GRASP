/*  Allocate memory and check for allocation failure; part of detect-lines.
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
#include <stdlib.h>


void *xmalloc(size)
  size_t size;
{
  void *ptr;

  if (size == 0)
    size = 1;
  ptr = malloc(size);
  if (ptr == NULL)
    print_error(ERR_NOMEM,NULL);
  return ptr;
}



void *xcalloc(nelem, elsize)
  size_t nelem;
  size_t elsize;
{
  void *ptr;

  if (elsize == 0)
    elsize = 1;
  if (nelem == 0)
    nelem = 1;
  ptr = calloc(nelem,elsize);
  if (ptr == NULL)
    print_error(ERR_NOMEM,NULL);
  return ptr;
}



void *xrealloc(old_ptr, size)
  void   *old_ptr;
  size_t size;
{
  void *ptr;

  if (size == 0)
    size = 1;
  ptr = realloc(old_ptr,size);
  if (ptr == NULL)
    print_error(ERR_NOMEM,NULL);
  return ptr;
}
