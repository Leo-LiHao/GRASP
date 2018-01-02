/* Test program for the regular expression routines */

#ifndef lint
static char *sccsid="@(#)retest.c	1.1 (LBL) 3/29/85";
static char sccsid[] = "@(#)retest.c	4.1	(Berkeley)	10/19/82";
static char Version[] =
   "$Id: retest.c,v 1.1.1.1 1997/03/18 00:44:37 cvs Exp $";
#endif

#include <stdio.h>
#include <ctype.h>

int	l_onecase = 0;
char	*_start;
char	*_escaped;
char	*convexp();
char	*expmatch();



main()
{
   char reg[132];
   char *ireg;
   char str[132];
   char *match;
   char matstr[132];
   char c;
   
   for (;;) {
      printf ("\nexpr: ");
      scanf ("%s", reg);
      ireg = convexp(reg);
      match = ireg;
      while(*match) {
	 switch (*match) {
	    
	 case '\\':
	 case '(':
	 case ')':
	 case '|':
	    printf ("%c", *match);
	    break;
	    
	 default:
	    if (isalnum(*match))
	       printf("%c", *match);
	    else
	       printf ("<%03o>", *match);
	    break;
	 }
	 match++;
      }
      printf("\n");
      getchar();
      for (;;) {
	 printf ("string: ");
	 match = str;
	 while ((c = getchar()) != '\n')
	    *match++ = c;
	 *match = 0;
	 if (str[0] == '#')
	    break;
	 matstr[0] = 0;
	 _start = str;
	 _escaped = 0;
	 match = expmatch (str, ireg, matstr);
	 if (match == 0)
	    printf ("FAILED\n");
	 else
	    printf ("match\nmatstr = %s\n", matstr);
      }
      
   }
}