/*  Print an error message and terminate the program; part of detect-lines.
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



/* Print an error message and terminate the program */
void print_error(message, info)
  const char *message;
  const char *info;
{
  if (info != NULL)
    fprintf(stderr,"%s: %s %s\n","detlines",message,info);
  else
    fprintf(stderr,"%s: %s\n","detlines",message);

  exit(1);
}
