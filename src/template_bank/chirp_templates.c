/* GRASP: Copyright 1997,1998  Bruce Allen */
static char *rcsid="$Id: chirp_templates.c,v 1.6 1998/02/04 01:33:12 ballen Exp $\
\n$Name: RELEASE_1_9_8 $";

#include "grasp.h"

 struct chirp_space global_space;

/* This global parameter is assigned by the routines
   get_chirp_templates(), to be used by the function chirp_metric().
   This allows chirp_metric() to access these parameters while still
   having a generic argument list.  This allows it to be passed to
   routines such as tiling_2d(), which doesn't need to know the
   details of how the parameter space metric is computed. */


void set_chirp_space(struct chirp_space space)

     /* This routine simply sets the global parameter global_space,
	above, so that the routine chirp_metric (below) can be used.

	The argument is:

	space: Input.  The structure that global_space is set to
	  equal. */

{
  global_space=space;
  return;
}


int chirp_metric(double x, double y, double *m)

     /* This routine computes the coefficients of a local distance
	metric on the space of chirp templates, at a specified
	location in that space.  It returns 0 upon successful
	completion, or 1 if the metric could not be computed at that
	point.

	The arguments are:

	x: Input.  The x coordinate of the specified point.

	y: Input.  The y coordinate of the specified point.

	m: Output.  The array m[0..2] contains the three independent
	  components of the local distance metric; see below.

	Comments: This routine uses the global parameter global_space,
	  defined at the top of this module.  This parameter must be
	  set before chirp_metric() can be called.

	The x,y coordinate system used is the one defined in the
	global_space structure: it is to a (counterclockwise) rotation
	of the tau0,tau1 coordinate system by an angle of
	global_space.angle.  The metric is computed by interpolating
	the grid of precomputed quadratic coefficients of the match
	function, stored in global_space.grid.

	The local distance function is defined so that the match
	decreases to the value global_space.match at a proper distance
	of 1.  That is, the unit of proper distance is defined to be
	the maximum template patch radius.  The metric coefficients
	m[0..2] are related to the proper interval ds^2 by:

	  ds^2 = m[0]*dx*dx + m[1]*dx*dy + m[2]*dy*dy

	Note that ds^2=1 corresponds to the match being equal to
	m=global_space.match in the definition of struct cubic_grid.
	We ignore the cubic components of the match function when
	dealing with the local distance metric.  So we have:

	  global_space.match = 1 + coef[0]*dx*dx +coef[1]*dx*dy
	                         + coef[2]*dy*dy

	when ds^2=1.  So the metric components m[] are related to the
	match function coefficients coef[] via:

	  m[i] = coef[i] / (global_space.match - 1),  i=0..2.

	This routine also checks to make sure that the resulting
	metric is positive definite, which should always be the case
	if the match function is locally paraboloidal (rather than a
	saddle). */

{
  double tau0,tau1,m_tot,eta,m1,m2;
  float coef[10];

  /* Convert from (x,y) to (m1,m2) coordinates. */
  tau0=x*cos(global_space.angle)-y*sin(global_space.angle);
  tau1=x*sin(global_space.angle)+y*cos(global_space.angle);
  if(!m_and_eta(tau0,tau1,&m_tot,&eta,global_space.m_mn,
		global_space.m_mx,M_PI*global_space.ftau))
    return 1;
  eta=pow(1.0-4.0*eta,0.5);
  m1=0.5*m_tot*(1.0-eta);
  m2=0.5*m_tot*(1.0+eta);

  /* Get the coefficients of the match function at (m1,m2). */
  if(get_cubic(m1,m2,global_space.grid,coef))
    return 1;

  /* Check that the metric will be positive definite. */
  if(coef[0]>0.0 || coef[2]>0.0 ||
     coef[1]*coef[1]>4.0*coef[0]*coef[2])
    return 1;

  /* Assign the output variables, and return. */
  m[0]=coef[0]/(global_space.match-1.0);
  m[1]=coef[1]/(global_space.match-1.0);
  m[2]=coef[2]/(global_space.match-1.0);
  return 0;
}


void get_chirp_boundary(struct chirp_space *space)

     /* This routine computes the boundary of the triangular parameter
	space of chirp signals, setting the fields x_bound and y_bound
	of *space.  Memory for these arrays is allocated in this
	routine; to free it, call free((*space).x_bound) and
	free((*space).y_bound).

	The argument is:

	space: Input/Output.  The data structure describing the
	  parameter space of chirps.  This routine uses as input the
	  fields m_mn, m_mx, ftau, angle, and n_bound, and assigns the
	  fields x_bound and y_bound.

	Comments: The boundary is actually made to lie just inside of
	  the triangular region specified by m_mn<m2<m1<m_mx.
	  Particular care is taken along the equal mass line to insure
	  that when one connects the computed points, the resulting
	  line segments lie entirely inside the true boundary: the
	  points must be shifted inwards to account for the slight
	  concavity along this curve.  Such caution is required
	  because we must occasionally convert points back into mass
	  space, and bad things happen if the points lie even slightly
	  past the equal mass line. */

{
  double m,dm,dm_long,*x,*y,*x_eq,*y_eq,x_tmp,pf,cos_angle,sin_angle;
  int i,n_side,n_longside;

  /* Compute the number and interval of points along each of the three
     sides of the space.  If (*space).n_bound is not divisible by 3,
     the sides will have different numbers of points; we place the
     extra ones on the lower mass line. */
  n_side=(*space).n_bound/3;
  n_longside=(*space).n_bound-2*n_side;
  dm=((*space).m_mx-(*space).m_mn)/n_side;
  dm_long=((*space).m_mx-(*space).m_mn)/n_longside;

  /* Compute pi times the frequency used to define tau0,tau1. */
  pf=M_PI*(*space).ftau;

  /* Allocate the output arrays. */
  x=(double *)malloc(((*space).n_bound+1)*sizeof(double));
  y=(double *)malloc(((*space).n_bound+1)*sizeof(double));
  (*space).x_bound=x;
  (*space).y_bound=y;

  /* Allocate temporary arrays for the equal mass line, excluding the
     corners. */
  x_eq=(double *)malloc((n_side-1)*sizeof(double))-1;
  y_eq=(double *)malloc((n_side-1)*sizeof(double))-1;

  /* Compute points along the equal mass line. */
  for(i=0,m=(*space).m_mn;i<n_side;i++,m+=dm)
    tau_of_mass(m,m,pf,x+i,y+i);

  /* Compute points along the upper mass line. */
  for(i=n_side,m=(*space).m_mx;i<2*n_side;i++,m-=dm)
    tau_of_mass((*space).m_mx,m,pf,x+i,y+i);

  /* Compute remaining points along the lower mass line. */
  for(i=2*n_side,m=(*space).m_mx;i<(*space).n_bound;i++,m-=dm_long)
    tau_of_mass(m,(*space).m_mn,pf,x+i,y+i);

  /* In tau space, the equal mass line is concave -- adjust the points
     inward so that the connecting line segments will inscribe the
     true curve. */
  for(i=1;i<n_side;i++){
    x_eq[i]=1.5*x[i]-0.25*(x[i-1]+x[i+1]);
    y_eq[i]=1.5*y[i]-0.25*(y[i-1]+y[i+1]);
  }
  for(i=1;i<n_side;i++){
    x[i]=x_eq[i];
    y[i]=y_eq[i];
  }

  /* The corner adjustments are potentially quite complicated.  For
     now, just round them off. */
  x[0]=0.5*(x[1]+x[(*space).n_bound-1]);
  y[0]=0.5*(y[1]+y[(*space).n_bound-1]);
  x[n_side]=0.5*(x[n_side+1]+x[n_side-1]);
  y[n_side]=0.5*(y[n_side+1]+y[n_side-1]);

  /* Now transform from tau0,tau1 coordinates to x,y coordinates. */
  cos_angle=cos((*space).angle);
  sin_angle=sin((*space).angle);
  for(i=0;i<(*space).n_bound;i++){
    x_tmp=x[i]*cos_angle+y[i]*sin_angle;
    y[i]=-x[i]*sin_angle+y[i]*cos_angle;
    x[i]=x_tmp;
  }

  /* Explicitly close off the curve. */
  x[(*space).n_bound]=x[0];
  y[(*space).n_bound]=y[0];

  /* Free the temporary arrays and return. */
  free(x_eq+1);
  free(y_eq+1);
  return;
}


int get_chirp_grid(struct chirp_space *space, const char *gridfile)

     /* This routine sets the field (*space).grid, which contains
	pre-computed coefficients of the cubic fit to the match
	function at various points over the parameter space.  It
	returns 0 if no warnings were generated, 1 if parameters used
	to generate the coefficient grid differed in some nontrivial
	way from those of the parameter space, or 2 if the coefficient
	grid could not be read in; in the latter case, (*space).grid
	is unchanged.  Otherwise, this routine will allocate memory
	for the coefficient grid; to free this memory, call
	free_cubic((*space).grid).

	The arguments are:

	space: Input/Output.  The parameter space over which the
	  coefficient grid is being assigned.  The fields m_mn and
	  m_mx are used only to check that the grid covers the space.
	  The fields ftau, angle, and match are used to check, rotate,
	  and rescale the grid's coordinate system (repsectively).
	  The field grid is the one which is set by this routine.

	gridfile: Input.  The name of a file containing the
	  pre-computed coefficients of the match function; see the
	  routines generate_cubic() and read_cubic(). */

{
  struct cubic_grid grid;
  int code=0;

  /* Read in the grid of match coefficients. */
  if(read_cubic(&grid,gridfile)){
    GR_start_error("get_chirp_grid()",rcsid,__FILE__,__LINE__);
    GR_report_error("Grid file %s is corrupt.\n",gridfile);
    GR_end_error();
    return 2;
  }

  /* Check whether the grid covers the parameter space. */
  if((*space).m_mn<grid.m_mn || (*space).m_mx>grid.m_mx){
    GR_start_error("get_chirp_grid()",rcsid,__FILE__,__LINE__);
    GR_report_error("Warning: Parameter space extends outside of the"
		    " grid of precomputed match\n"
		    "coefficients.\n"
		    "Parameter space mass range:  [%f,%f]\n"
		    "Coefficient grid mass range: [%f,%f]\n",
		    (*space).m_mn,(*space).m_mx,grid.m_mn,grid.m_mx);
    GR_end_error();
    code=1;
  }

  /* Check whether the grid was computed with the same starting
     frequency for the mass-to-tau conversion.  If not, the one in the
     grid data structure will be changed; hopefully the match
     coefficients are not too sensitively dependent on it. */
  if((*space).ftau!=grid.ftau){
    GR_start_error("get_chirp_grid()",rcsid,__FILE__,__LINE__);
    GR_report_error("Warning: Parameter space starting frequency does"
		    " not agree with that used to\n"
		    "compute the match coefficients.\n"
		    "Parameter space  ftau=%f\n"
		    "Coefficient grid ftau=%f\n",(*space).ftau,
		    grid.ftau);
    GR_end_error();
    grid.ftau=(*space).ftau;
    code=1;
  }

  /* Rotate and rescale the match coefficients so that they agree with
     the current parameter space's coordinates and match level. */
  transform_cubic(&grid,(*space).angle,(*space).match);

  /* Assign the (*space).grid field and return the error code. */
  (*space).grid=grid;
  return code;
}


int get_chirp_templates(struct chirp_space *space)

     /* This routine computes the positions of a mesh of chirp
	templates on a parameter space, using the generic tiling
	routine tiling_2d().  It passes to tiling_2d() the boundary
	polygon defined by (*space).x_bound, (*space).y_bound, and
	(*space).n_bound, as well as the template space metric
	function chirp_metric(), then converts the returned list of
	templates into the array (*space).templates.  The
	chirp_metric() routine requires the parameter space to be
	passed to it as a global static variable named global_space;
	this variable is set equal to *space.  The routine
	get_chirp_templates() itself returns the error code generated
	by the call to tiling_2d(); see the documentation of that
	routine.  If a fatal error occurs before tiling_2d() is even
	called, get_chirp_templates() returns an error code of 3,
	(*space).n_templates is set to 0 and (*space).templates to
	NULL.  Otherwise, memory for the template array will be
	allocated; to free this memory, call free((*space).templates).

	The argument is:

	space: Input/Output.  The data structure of the parameter
	  space being filled with templates.  This routine uses as
	  input the fields n_bound, x_bound, y_bound, and grid, and
	  assigns the fields n_templates and templates. */

{
  struct tile *head,*here;
  struct chirp_template *t;
  double tmp,cos_angle,sin_angle;
  int i,code;

  /* Set default return values, and assign the global variable. */
  (*space).templates=NULL;
  (*space).n_templates=0;
  set_chirp_space(*space);

  /* Make sure that the parameter space boundary has been defined. */
  if((*space).x_bound==NULL || (*space).y_bound==NULL ||
     (*space).n_bound<=0){
    GR_start_error("get_chirp_templates()",rcsid,__FILE__,__LINE__);
    GR_report_error("Parameter space boundary has not been set;\n"
		    "see the routine get_chirp_boundary().\n");
    GR_end_error();
    return 3;
  }

  /* Make sure that the grid of match coefficients has been set, and
     assign this to the global variable grid. */
  if((*space).grid.coef==NULL || (*space).grid.n==0){
    GR_start_error("get_chirp_templates()",rcsid,__FILE__,__LINE__);
    GR_report_error("Match coefficient grid has not been set;\n"
		    "see the routine get_chirp_grid().\n");
    GR_end_error();
    return 3;
  }

  /* Now allocate a dummy handle to template linked list. */
  head=here=(struct tile *)malloc(sizeof(struct tile));

  /* Generate the linked list of templates. */
  code=tiling_2d((*space).x_bound,(*space).y_bound,(*space).n_bound,
		 chirp_metric,&here,&(*space).n_templates);

  /* Free the dummy handle. */
  here=head;
  head=(*here).next;
  free(here);

  /* Compute some temporary variables. */
  cos_angle=cos((*space).angle);
  sin_angle=sin((*space).angle);

  /* Convert the list into a chirp_template array, freeing the list as
     we go.  If tiling_2d() is working properly, no errors should
     occur, but check anyhow. */
  (*space).templates=(struct chirp_template *)
    malloc((*space).n_templates*sizeof(struct chirp_template));
  for(i=0,here=head;i<(*space).n_templates && here!=NULL;
      i++,here=head){

    /* Simplify notation with a temporary variable. */
    t=(*space).templates+i;

    /* Assign template fields which have already been calculated by
       tiling_2d(). */
    (*t).num=i;
    (*t).flag=(*here).flag;
    (*t).x=(*here).x;
    (*t).y=(*here).y;
    (*t).dx=(*here).dx;
    (*t).dy=(*here).dy;
    (*t).semimajor=(*here).r1;
    (*t).semiminor=(*here).r2;
    (*t).theta=(*here).theta;

    /* Compute template location in tau0,tau1 coordinates. */
    (*t).tau0=(*t).x*cos_angle-(*t).y*sin_angle;
    (*t).tau1=(*t).x*sin_angle+(*t).y*cos_angle;

    /* Compute template location in various mass coordinates. */
    if(m_and_eta((*t).tau0,(*t).tau1,&((*t).mtotal),&((*t).eta),
		 (*space).m_mn,(*space).m_mx,M_PI*(*space).ftau)){
      (*t).mred=(*t).eta*(*t).mtotal;
      (*t).mchirp=(*t).mred*pow((*t).eta,-0.4);
      tmp=pow(1.0-4.0*(*t).eta,0.5);
      (*t).m1=0.5*(*t).mtotal*(1.0-tmp);
      (*t).m2=0.5*(*t).mtotal*(1.0+tmp);
    }else{
      GR_start_error("get_chirp_templates()",rcsid,__FILE__,__LINE__);
      GR_report_error("Template %i is out of bounds, at"
		      " (tau0,tau1)=(%f,%f)\n",i,(*t).tau0,(*t).tau1);
      GR_end_error();
      if(!((*t).flag))
	(*t).flag=-1;
    }

    /* Free array as we go. */
    head=(*here).next;
    free(here);
  }

  /* Check whether the list really was of length n_templates. */
  for(here=head;here!=NULL;here=head,i++){
    head=(*here).next;
    free(here);
  }
  if(i!=(*space).n_templates){
    GR_start_error("get_chirp_templates()",rcsid,__FILE__,__LINE__);
    GR_report_error("List contained %i tiles instead of %i.\n",i,
		    (*space).n_templates);
    GR_end_error();
  }

  /* Return the error code generated by template_2d(). */
  return code;
}


int plot_chirp_templates(struct chirp_space space,
			 double magnification, int plot_boundary,
			 int plot_patches, int plot_ellipses,
			 int plot_flags, const char *psfile)

     /* This routine plots a postscript file of a chirp parameter
	space and the templates covering it, using the generic
	template-plotting routine plot_list().  The input parameters
	for plot_list() are passed directly from the parameters for
	plot_chirp_templates(), except as follows: the boundary is
	taken from space.x_bound, space.y_bound, and space.n_bound,
	the linked list of templates is generated from space.templates
	and space.n_templates, and the rotation angle of the plot is
	such that tau0 runs down the length of the page.  The routine
	plot_chirp_templates() itself returns the number of pages of
	postscript output.

	The arguments are:

	space: Input.  The data structure of the parameter space and
	  its templates.  This routine makes use of the fields
	  n_bound, x_bound, y_bound, and angle; if template patches or
	  flags are to be plotted, the fields n_templates and
	  templates are also used.

	magnification: Input.  The scale factor of points (1/72 of an
	  inch) per unit coordinate distance in the parameter space.

	plot_boundary: Input.  1 if boundary is to be shown, 0
	  otherwise.

	plot_patches: Input.  1 if rectangular brickwork of patches is
	  to be shown, 0 otherwise.

	plot_ellipses: Input.  1 if the overlapping match contours are
	  to be shown, 0 otherwise.

	plot_flags: Input.  If nonzero, indicates the size of dot used
	  to mark flagged templates (in points = 1/72 inches).  If
	  zero, flags are ignored.

	psfile: Input.  The name of the postscript file created. */

{
  struct tile *here,*head;
  int i,code;

  /* If necessary, generate a linked list of tiles from the template
     array. */
  here=head=(struct tile *)malloc(sizeof(struct tile));
  if(plot_patches || plot_ellipses || plot_flags)
    for(i=0;i<space.n_templates;i++){
      (*here).next=(struct tile *)malloc(sizeof(struct tile));
      here=(*here).next;
      (*here).flag=space.templates[i].flag;
      (*here).x=space.templates[i].x;
      (*here).y=space.templates[i].y;
      (*here).dx=space.templates[i].dx;
      (*here).dy=space.templates[i].dy;
      (*here).r1=space.templates[i].semimajor;
      (*here).r2=space.templates[i].semiminor;
      (*here).theta=space.templates[i].theta;
    }
  (*here).next=NULL;

  /* Generate postscript using plot_list(). */
  code=plot_list(space.x_bound,space.y_bound,space.n_bound,
		 (*head).next,space.n_templates,space.angle-0.5*M_PI,
		 magnification,plot_boundary,plot_patches,
		 plot_ellipses,plot_flags,psfile);

  /* Free the linked list. */
  for(here=head;here!=NULL;here=head){
    head=(*here).next;
    free(here);
  }

  /* Return the number of pages of postscript output. */
  return code;
}
