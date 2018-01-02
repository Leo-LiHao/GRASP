/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"

int main(int argc, char **argv)
{
  struct chirp_space space;
  double mag;
  int i,code,plot_boundary,plot_patches,plot_ellipses,plot_flags;
  FILE *fpout;

  /* Set the chirp space parameters. */
  space.m_mn=1.0;
  space.m_mx=3.0;
  space.angle=0.318;
  space.ftau=140.0;
  space.match=0.98;
  space.n_bound=600;

  /* Set the plotting parameters. */
  mag=4200.0;         /* Magnification factor. */
  plot_boundary=1;    /* Do show the boundary. */
  plot_patches=1;     /* Do show the patches. */
  plot_ellipses=0;    /* Don't show the circumscribed ellipses. */
  plot_flags=1;       /* Do indicate any flagged templates. */

  /* Generate the parameter space boundary. */
  get_chirp_boundary(&space);

  /* Print out the boundary points. */
  fpout=fopen("boundary.ascii","w");
  for(i=0;i<=space.n_bound;i++)
    fprintf(fpout,"%f %f\n",space.x_bound[i],space.y_bound[i]);
  fclose(fpout);

  /* Get the match function parameters over the space. */
  if(get_chirp_grid(&space,"cubic_coef_40meter_m=0.8-3.2.ascii"))
    return 1;

  /* Generate the template mesh. */
  code=get_chirp_templates(&space);
  fprintf(stdout,"%s: get_chirp_templates generated %i templates and"
	  " exited with code %i.\n",argv[0],space.n_templates,code);

  /* Print a list of template positions in mass space. */
  fpout=fopen("templates.ascii","w");
  for(i=0;i<space.n_templates;i++)
    fprintf(fpout,"%f %f\n",space.templates[i].m1,
	    space.templates[i].m2);
  fclose(fpout);

  /* Flag every template, so that their centres will be marked with
     dots on the plot.  Normally, a template will be flagged only if
     an error occured while generating it. */
  for(i=0;i<space.n_templates;i++)
    space.templates[i].flag=1;

  /* Plot a postscript diagram of the parameter space. */
  code=plot_chirp_templates(space,mag,plot_boundary,plot_patches,
			    plot_ellipses,plot_flags,"templates.ps");
  fprintf(stdout,"%s: plot_chirp_templates generated %i pages of"
	  " postscript output.\n",argv[0],code);

  return 0;
}
