/* GRASP: Copyright 1997,1998  Bruce Allen */
static char *rcsid="$Id: tiling_2d.c,v 1.7 1998/10/05 22:46:41 ballen Exp $\
\n$Name: RELEASE_1_9_8 $";

#include "grasp.h"

#define RtD 57.2957795131 /* Conversion factor radians -> degrees */
#define N_COLS 1000000
#define N_ROWS 1000000
#define X_MARGIN 36
#define Y_MARGIN 36
#define X_SIZE 612
#define Y_SIZE 792
#define N_OBJ 797

int tiling_2d(double *x_bound, double *y_bound, int npts,
	      int (*metric)(double , double , double *),
	      struct tile **tail, int *n_tiles)

     /* This is a generic routine for laying out a mesh of overlapping
        rectangular tiles in a (relatively) arbitrary two-dimensional
        parameter space.  The tiles are sized such that no point on
        the tile is more than one unit of proper distance from its
        centre, where the proper distance is computed with a metric
        function which can vary arbitrarily over the parameter space.
        Thus the size and shape of the tiles can and will vary over
        the space.  The tiles are rectangular, aligned with the
        coordinate axes, and are laid out in a simple brickwork
        lattice, with extra overlapping tiles on the edges to ensure
        complete coverage of the space.  The lattice of tile positions
        is stored as a linked list.

	Note that the routine can be easily modified to lay out tiles
	with non-unit radius.  To scale the tiles by a factor d,
	simply multiply the metric function by a factor 1/(d*d).

	Some details of the algorithm: The routine operates by
	dividing the parameter space into columns, and then calling
	the routine column() to fill each column.  The column widths
	are calculated so as to maximize the individual tile areas
	within that column, but the height of each tile may vary
	within each column, and the width may between one column and
	the next.  The routine corner() places extra tiles at the
	corners of each column, so as to completely cover the boundary
	of the space.  The result is a lattice of tile positions such
	that (1) no point in the space is more than unit proper
	distance from the centre of a tile, and (2) no tile is centred
	on a point outside the boundary of the space.  Within these
	restrictions, every effort is made to minimize the number of
	tiles, by maximizing the proper area of each tile.

	Upon successful execution, tiling_2d() attaches the new linked
	list to (**tail).next, updates *tail to point to the new tail
	of the list, and returns a value of 0.  It returns an error
	code of 1 if it suspects that some of the columns may not be
	properly filled.  It returns 2 if at any point the width of
	the parameter space was more than N_COLS times the computed
	column width.  It returns 3 if the routine terminated
	prematurely for any other reason (usually because the
	algorithm accidentally stepped out of the parameter space, due
	to imprecise interpolation of the boundary).  In the case of
	error codes 2 and 3, tiling_2d() still attaches the generated
	list onto (**tail).next (up to the point where the error
	occurred), but does not update the position of *tail.
	tiling_2d() will also write appropriate error messages
	indicating where any errors took place, and the flag field of
	the tile on the list (at the time of the error) is set to the
	error code.  A flag value of -1 on any tile is a warning
	flag; the tile was placed correctly, but not all of the fields
	were calculated.

	The arguments are:

	x_bound: Input.  An array [0..npts] storing the x-coordinates
	  of npts points along the boundary of the parameter space.
	  The array has length npts+1, but the index [npts] should
	  refer to the same point as [0].

	y_bound: Input.  As above, but the y-coordinates.

	npts: Input.  The number of points used to specify the
	  boundary.

	metric(): Input.  This function computes the three independent
	  components of the distance metric matrix at a given point in
	  parameter space.  The first two arguments are the x and y
	  coordinates of the requested point, the third passes back
	  the metric components in a three-element array.  The [0],
	  [1], and [2] metric components are defined in terms of the
	  proper interval as follows:

	    ds^2 = [0]*dx*dx + [1]*dx*dy + [2]*dy*dy .

	  metric() itself should return 0, or 1 if the metric is
	  undefined or not computable at the specified location.

	tail: Input/Output.  Initially points to the ``tail'' of a
	  pre-existing list; the generated list is attached to
	  (**tail).next.  Upon successful completion, **tail is
	  updated to the new tail of the list.

	n_tiles: Input/Output.  A running tally of the number of tiles
	  in the mesh; it is incremented each time a tile is added
	  (and decremented whenever a tile is removed).

	The tiling_2d() routine makes very few assumptions about the
	parameter space.  The most stringent is the assumption that
	the parameter space boundary can be expressed as bivalued
	functions of both x and y; that is, both vertical and
	horizontal lines intersect the boundary at no more than two
	points.  If a vertical line intersects at more than two
	points, the routine may come to a point where it cannot
	determine the location or width of the next column, and will
	terminate.  If a horizontal line intersects at more than two
	points, the routine column() may encounter difficulties
	filling up the ``corners'' of the columns, which may result in
	incomplete coverage of the parameter space.  These conditions
	are checked by the routine test_cycle(). */

{
  double left,colwidth,tempwidth,bottom,top,leftmost,rightmost;
  int code,return_code=0;
  struct tile *here,*temp1,*temp2;

  /* Make sure there is a pre-existing list to attach to. */
  if(*tail==NULL){
    GR_start_error("tiling_2d()",rcsid,__FILE__,__LINE__);
    GR_report_error("tail pointer must be allocated before being"
		    " passed.\n");
    GR_end_error();
    return 3;
  }
  here=*tail;

  /* Check to see whether the boundary is likely to give any problems.
     If it is, continue anyway, but print a warning. */
  code=test_cycle(x_bound,npts,&leftmost,&rightmost);
  if(code!=2){
    GR_start_error("tiling_2d()",rcsid,__FILE__,__LINE__);
    GR_report_error("Vertical lines can intersect the boundary at"
		    " more than two points.\n"
		    "This may cause difficulties in locating the top"
		    " and bottom of a column.\n");
    GR_end_error();
  }
  code=test_cycle(x_bound,npts,&leftmost,&rightmost);
  if(code!=2){
    GR_start_error("tiling_2d()",rcsid,__FILE__,__LINE__);
    GR_report_error("Horizontal lines can intersect the boundary at"
		    " more than two points.\n"
		    "This may cause difficulties in placing tiles"
		    " near the corners of a column.\n");
    GR_end_error();
  }

  /* Set up the position of the column ``previous'' to the first one
     (since the main loop assumes it is adding a new column after a
     previous one). */
  left=leftmost;
  colwidth=0.0;
  code=0;

  while((left+=colwidth)<rightmost){

    /* Compute the column width.  If the previous attempt at placing
       the column returned a code of 2, then most likely the column
       was too wide; reduce colwidth and try again. */
    if(code==2)
      colwidth*=0.5;

    /* Otherwise, set the column width to be the smaller of the widths
       computed at the lower left and upper left corners. */
    else{
      if(range(x_bound,y_bound,npts,left,&bottom,&top)!=2){
	GR_start_error("tiling_2d()",rcsid,__FILE__,__LINE__);
	GR_report_error("Could not compute column height at x=%f\n",
			left);
	GR_end_error();
	return (*here).flag=3;
      }
      if(get_width(metric,left,bottom,&colwidth)){
	GR_start_error("tiling_2d()",rcsid,__FILE__,__LINE__);
	GR_report_error("Could not compute column width at"
			" (x,y)=(%f,%f)\n",left,bottom);
	GR_end_error();
	return (*here).flag=3;
      }
      if(get_width(metric,left,top,&tempwidth)){
	GR_start_error("tiling_2d()",rcsid,__FILE__,__LINE__);
	GR_report_error("Could not compute column width at"
			" (x,y)=(%f,%f)\n",left,top);
	GR_end_error();
	return (*here).flag=3;
      }
      if(tempwidth<colwidth)
	colwidth=tempwidth;
    }

    /* If the computed column width is too small, return an error of
       type 2. */
    if(N_COLS*colwidth<rightmost-leftmost){
      GR_start_error("tiling_2d()",rcsid,__FILE__,__LINE__);
      GR_report_error("Projected number of columns exceeded %i at"
		      " x=%f\n",N_COLS,left);
      GR_end_error();
      return (*here).flag=2;
    }

    /* Make sure that the column doesn't hang off the rightmost end of
       the parameter space, by reducing the column width if necessary.
       This last column is permitted to be narrower than the parameter
       space width / N_COLS, since the next-to-last column might just
       happen to end very near the rightmost point. */
    if(left+colwidth>rightmost)
      colwidth=rightmost-left;

    /* Now lay down the column. */
    code=column(x_bound,y_bound,npts,left,colwidth,metric,&here,
		n_tiles);

    /* A return code of 1 from column() means that the routine was
       able to place at least some tiles in the column, but encoutered
       difficulties in placing tiles near the boundaries.  Normally
       this is because various numerical approximations caused it to
       overstep the boundaries.  This does not terminate the routine,
       but it may result in incomplete coverage around the
       boundaries. */
    if(code==1)
      return_code=1;

    /* If the heights of tiles in a column ever dropped below the
       column height / N_ROWS, column() returns an error code 2.
       Normally this is because the requested column width was to
       large, so de-allocate the failed column and try again.  (The
       column width will be reduced at the start of the loop, where
       colwidth is computed.) */
    else if(code==2){
      for(temp1=temp2=(*here).next;temp1!=NULL;temp2=temp1){
	temp1=(*temp1).next;
	free(temp2);
	(*n_tiles)--;
      }
      (*here).next=NULL;
      /* The top of the loop automatically increments the value of
         left, so we'd better decrement it here so that it stays the
         same in the next loop. */
      left-=colwidth;
    }

    /* An error code of 3 from column() is a terminating error -- it
       could not even start the column successfully.  A decsriptive
       error message will already have been printed. */
    else if(code==3){
      GR_start_error("tiling_2d()",rcsid,__FILE__,__LINE__);
      GR_report_error("Received return code 3 from column()."
		      "  Exiting.\n");
      GR_end_error();
      return 3;
    }
  }

  /* Normal termination.  Point *tail to the new list tail, and return
     a code of 0 (flawless execution) or 1 (possible incomplete
     coverage). */
  *tail=here;
  return return_code;
}


int column(double *x_bound, double *y_bound, int npts, double left,
	   double width, int (*metric)(double , double , double *),
	   struct tile **tail, int *n_tiles)

     /* This routine places tiles in a column which completely covers
	the region in parameter space between the two x values left
	and left+width.  The new tiles are attached to (**tail).next,
	and *tail is set to point to the last tile in the list.  A
	return code of 0 indicates successful execution.  A return
	code of 1 means that the routine encountered difficulties
	placing tiles near the boundaries (possibly because numerical
	imprecision caused it to overstep the boundaries), but it was
	able to fill the column at least partially.  A return code of
	2 means that at some point the height of the tiles dropped to
	less than the column height / N_ROWS, probably because the
	requested column width was too large; in this case *tail is
	left unchanged.  A return code of 3 indicates a more serious
	error preventing column() from even starting the column;
	again, *tail is left unchanged.  Whenever an error of type 1
	or 3 occurs, the .flag field of the last tile on the list is
	set to the error code.  The .flag field is set to -1 on any
	``incomplete'' tile, for which the .r1, .r2, and .theta fields
	could not be computed.

	The arguments are:

	x_bound: Input.  An array [0..npts] storing the x-coordinates
	  of npts points along the boundary of the parameter space.
	  The array has length npts+1, but the index [npts] should
	  refer to the same point as [0].

	y_bound: Input.  As above, but the y-coordinates.

	npts: Input.  The number of points used to specify the
	  boundary.

	left: Input.  The x-value of the left side of the column.

	width: Input.  The width of the column.

	metric(): Input.  This function computes the three independent
	  components of the distance metric matrix at a given point in
	  parameter space.  See the documentation for tiling_2d(),
	  above.

	tail: Input/Output.  Initially points to the last tile of the
	  previous column.  Upon successful completion, **tail is
	  updated to the last tile of the current column.

	n_tiles: Input/Output.  A running tally of the number of tiles
	  in the mesh; it is incremented each time a tile is added.

	In addition to simply stacking tiles up the length of the
	column, column() calls the routine corner() to cover the
	``corners'' of the column which might otherwise be missed by a
	straightforward stacking.  These tiles may overlap
	significantly with tiles in this and adjacent columns, but at
	least they are guaranteed to cover the entire parameter space
	between the specified values of x.  A failure in corner() will
	result in an error of type 1. */

{
  double centre,bottom,top,y,dy;
  int code,return_code=0;
  struct tile *here,*last;

  /* Find the base and the top of the column. */
  centre=left+0.5*width;
  if(range(x_bound,y_bound,npts,centre,&bottom,&top)!=2){
    GR_start_error("column()",rcsid,__FILE__,__LINE__);
    GR_report_error("Could not compute column height at x=%f\n",
		    centre);
    GR_end_error();
    return (**tail).flag=3;
  }

  /* Place the column base. */
  code=get_height(metric,centre,bottom,width,&dy);
  if(code==1){
    GR_start_error("column()",rcsid,__FILE__,__LINE__);
    GR_report_error("Stepped out of bounds at (x,y)=(%f,%f)\n",centre,
		    bottom);
    GR_end_error();
    return (**tail).flag=3;
  }else if(N_ROWS*dy<top-bottom){
    GR_start_error("column()",rcsid,__FILE__,__LINE__);
    GR_report_error("Projected number of rows exceeded %i at"
		    " (x,y)=(%f,%f)\n",N_ROWS,centre,bottom);
    GR_end_error();
    return 2;
  }
  here=last=(**tail).next=(struct tile *)malloc(sizeof(struct tile));
  (*n_tiles)++;
  (*here).x=centre;
  (*here).y=bottom;
  (*here).dx=width;
  (*here).dy=dy;
  (*here).next=NULL;
  (*here).r1=(*here).r2=(*here).theta=0.0;
  if(((*here).flag=-get_ellipse(metric,centre,bottom,&((*here).r1),
				&((*here).r2),&((*here).theta)))){
    GR_start_error("column()",rcsid,__FILE__,__LINE__);
    GR_report_error("Warning: Ellipse contour undefined at"
		    " (x,y)=(%f,%f)\n",centre,bottom);
    GR_end_error();
  }

  /* Fill in the lower-left corner of the column. */
  if(corner(x_bound,y_bound,npts,-1,-1,metric,here,&last,n_tiles))
    return_code=1;

  /* Fill in the lower-right corner of the column. */
  if(corner(x_bound,y_bound,npts,1,-1,metric,here,&last,n_tiles))
    return_code=1;

  /* Fill up the length of the column. */
  while(top>(*here).y+0.5*(*here).dy){

    /* Compute the position of the next tile. */
    if(top<(y=(*here).y+(*here).dy))
      y=top;

    /* Compute tile height at that position. */
    code=get_height(metric,centre,y,width,&dy);
    if(code==1){
      GR_start_error("column()",rcsid,__FILE__,__LINE__);
      GR_report_error("Stepped out of bounds at (x,y)=(%f,%f)\n",
		      centre,y);
      GR_end_error();
      return (*here).flag=1;
    }else if(N_ROWS*dy<top-bottom){
      GR_start_error("column()",rcsid,__FILE__,__LINE__);
      GR_report_error("Projected number of rows exceeded %i at"
		      " (x,y)=(%f,%f)\n",N_ROWS,centre,y);
      GR_end_error();
      return 2;
    }

    /* Revise the position of the next tile.
    if(top<(y=(*here).y+0.5*(*here).dy+0.5*dy))
      y=top; */

    /* Allocate and attach the tile to the list. */
    here=(*last).next=(struct tile *)malloc(sizeof(struct tile));
    last=here;
    (*n_tiles)++;
    (*here).x=centre;
    (*here).y=y;
    (*here).dx=width;
    (*here).dy=dy;
    (*here).next=NULL;
    (*here).r1=(*here).r2=(*here).theta=0.0;
    if(((*here).flag=-get_ellipse(metric,centre,y,&((*here).r1),
				  &((*here).r2),&((*here).theta)))){
      GR_start_error("column()",rcsid,__FILE__,__LINE__);
      GR_report_error("Warning: Ellipse contour undefined at"
		      " (x,y)=(%f,%f)\n",centre,y);
      GR_end_error();
    }
  }

  /* Fill in the upper-left corner of the column. */
  if(corner(x_bound,y_bound,npts,-1,1,metric,here,&last,n_tiles))
    return_code=1;

  /* Fill in the upper-right corner of the column. */
  if(corner(x_bound,y_bound,npts,1,1,metric,here,&last,n_tiles))
    return_code=1;

  /* If everything worked, reposition *tail and return. */
  *tail=last;
  return return_code;
}


int corner(double *x_bound, double *y_bound, int npts, int xsign,
	   int ysign, int (*metric)(double , double , double *),
	   struct tile *here, struct tile **tail, int *n_tiles)

     /* This routine places overlapping tiles to fill in a corner of a
	column.  Depending on the values of xsign and ysign, the tiles
	will be placed off to the upper right, lower left, etc.  The
	tile *here defines the current top or bottom of the column
	(depending on ysign), and is the point from which the
	corner-branch of tiles is ``grown''; *here is also used as an
	internal variable, though of course its revised value is not
	passed back to the calling routine.  Instead, the branch is
	attached to (**tail).next, and upon completion *tail is set to
	point to the new tail of the list (the cornermost tile).
	corner() returns 0 if it executes successfully.  It returns 1
	if at any point it could not locate the next tile position,
	(usually because imprecise interpolation of the boundary
	caused it to step outside of the parameter space), or if the
	computed tile height was less than the column height / N_ROWS.
	These are not critical errors; the program will simply stop
	trying to fill in that corner.  However, this may result in
	incomplete coverage of the parameter space.  Whenever such an
	error occurs, the .flag field of the last tile on the list is
	set to 1.

	The arguments are:

	x_bound: Input.  An array [0..npts] storing the x-coordinates
	  of npts points along the boundary of the parameter space.
	  The array has length npts+1, but the index [npts] should
	  refer to the same point as [0].

	y_bound: Input.  As above, but the y-coordinates.

	npts: Input.  The number of points used to specify the
	  boundary.

	xsign: Input.  Either +/- 1, indicating whether the corner is
	  on the right/left side of the column.

	ysign: Input.  Either +/- 1, indicating whether the corner is
	  on the top/bottom of the column.

	metric(): Input.  This function computes the three independent
	  components of the distance metric matrix at a given point in
	  parameter space.  See the documentation for tiling_2d(),
	  above.

	here: Input.  The base tile at the top/bottom of the column
	  from which the corner branch will extend.  Also used
	  internally as a local variable, but changes to it are not
	  returned.

	tail: Input/Output.  The last tile of the list.

	n_tiles: Input/Output.  A running tally of the number of tiles
	  in the mesh; it is incremented each time a tile is added. */

{
  double x_corner,y_corner,x,y,temp,width,height,dy;
  int code;

  /* Find the extreme location of the corner. */
  width=(*here).dx;
  x_corner=(*here).x+0.5*xsign*width;
  if(range(x_bound,y_bound,npts,x_corner,&temp,&y_corner)!=2){
    GR_start_error("corner()",rcsid,__FILE__,__LINE__);
    GR_report_error("Could not compute column height at x=%f\n",
		    x_corner);
    GR_end_error();
    return (*here).flag=1;
  }
  height=y_corner-temp;
  if(ysign<0)
    y_corner=temp;

  /* As long as the corner lies above/below the current top/bottom of
     the column... */
  while(ysign*y_corner>ysign*(*here).y+0.5*(*here).dy){

    /* Check whether the corner lies above/below the height y of the
       centre of the next tile to be added. */
    if(ysign*y_corner>ysign*(y=(*here).y+ysign*(*here).dy)){

      /* If so, position the tile at the point where y intersects the
         boundary. */
      if(range(y_bound,x_bound,npts,y,&x,&temp)!=2){
	GR_start_error("corner()",rcsid,__FILE__,__LINE__);
	GR_report_error("Could not compute intersection of boundary"
			" with y=%f\n",y);
	GR_end_error();
	return (*here).flag=1;
      }
      if(xsign<0)
	x=temp;
    }

    /* If not, position the tile at the corner. */
    else{
      x=x_corner;
      y=y_corner;
    }

    /* Compute the tile height. */
    code=get_height(metric,x,y,width,&dy);
    if(code==1){
      GR_start_error("corner()",rcsid,__FILE__,__LINE__);
      GR_report_error("Stepped out of bounds at (x,y)=(%f,%f)\n",x,y);
      GR_end_error();
      return (*here).flag=1;
    }else if(N_ROWS*dy<height){
      GR_start_error("corner()",rcsid,__FILE__,__LINE__);
      GR_report_error("Projected number of rows exceeded %i at"
		      " (x,y)=(%f,%f)\n",N_ROWS,x,y);
      GR_end_error();
      return (*here).flag=1;
    }

    /* Allocate and attach the next tile to the list. */
    here=(**tail).next=(struct tile *)malloc(sizeof(struct tile));
    *tail=here;
    (*n_tiles)++;
    (*here).x=x;
    (*here).y=y;
    (*here).dx=width;
    (*here).dy=dy;
    (*here).next=NULL;
    (*here).r1=(*here).r2=(*here).theta=0.0;
    if(((*here).flag=-get_ellipse(metric,x,y,&((*here).r1),
				  &((*here).r2),&((*here).theta)))){
      GR_start_error("corner()",rcsid,__FILE__,__LINE__);
      GR_report_error("Warning: Ellipse contour undefined at"
		      " (x,y)=(%f,%f)\n",x,y);
      GR_end_error();
    }
  }
  return 0;
}


int test_cycle(double *cycle, int npts, double *low, double *high)

     /* This routine examines a closed sequence of points, and
	determines the maximum and minimum values in the sequence.  It
	returns of times that the sequence reverses direction.  This
	number is always even, and always positive unless the sequence
	is completely degenerate (a set of equal points).

	The arguments are:

	cycle: Input.  The array cycle[0..npts] contains the sequence
	  of points.  Note that the array is of length npts+1; the [0]
	  and [npts] index values refer to the same point, so the
	  array explicitly describes a closed cycle.  This convention
	  simplifies many algorithms, including this one.

	npts: Input.  The number of points in the cycle.

	low: Output.  The minimum value in the cycle.

	high: Output.  The maximum value in the cycle.

	n_reverse: Output.  The number of times that the cycle
	  reverses direction.

	Comments: In the context of the tiling_2d() routine, this can
	  be used on the arrays containing the x (and y) coordinates
	  of points along the boundary, to find the total domain (and
	  range) of the parameter space.  Also, the number of
	  direction reversals is, equivalently, the maximum number of
	  times that a vertical (or horizontal) line will intersect
	  the boundary.  If this number is more than 2, the tiling
	  algorithm may encounter difficulties.

	  A note on conventions: If the sequence contains a set of
	  equal elements in adjacent positions, it is treated as
	  ``increasing'' through this set (i.e. it is as if subsequent
	  elements were slightly larger than the earlier elements).
	  Thus a sequence that is generally increasing with occasional
	  stops is treated as always increasing, whereas a sequence
	  that is decreasing with occasional stops is treated as
	  reversing direction twice at each stop. */

{
  int i,increasing,n_reverse;

  /* Set some preliminary values. */
  *low=*high=cycle[0];
  n_reverse=0;
  increasing=(cycle[npts]>=cycle[npts-1]);

  /* Loop through the cycle, updating the values of the output
     variables as we go. */
  for(i=0;i<npts;i++){
    if(cycle[i]<*low)
      *low=cycle[i];
    else if(cycle[i]>*high)
      *high=cycle[i];
    if(increasing!=(cycle[i+1]>=cycle[i])){
      increasing=!increasing;
      n_reverse++;
    }
  }
  return n_reverse;
}


int range(double *x_bound, double *y_bound, int npts, double x,
	  double *y_low, double *y_high)

     /* This function computes the intersection points between a
	closed boundary curve and a specified vertical line, by
	interpolating between points on the boundary.  It returns the
	number of intersection points found.

	If no intersection points are found, then y_low and y_high are
	left unchanged.  If one is found, they are both set to that
	value (this happens, for instance, if the vertical line just
	touches the edge of the boundary).  Otherwise, they are set to
	the lowest and highest points of intersection.

	The arguments are:

	x_bound: Input.  The array x_bound[0..npts] contains the x
	  components of a set of npts boundary points.  Note that the
	  array is of length npts+1; the [0] and [npts] index values
	  refer to the same point, so the array explicitly describes a
	  closed boundary.  This convention simplifies many
	  algorithms.

	y_bound: Input.  The array y_bound[0..npts] contains the y
	  components of the boundary points, as above.

	npts: Input.  The number of points along the boundary.

	x: Input.  The position of the vertical line.

	y_low: Output.  The y-value of the lowest intersection point.

	y_high: Output.  The y-value of the highest intersection
	  point.

	Comments: By interchanging x's and y's, this same routine will
	  find the intersections between the boundary and a horizontal
	  line.  In fact, the tiling_2d() algorithm calls the routine
	  more frequently in this manner. */

{
  double y;
  int i,n_int=0,increasing;

  /* The variable increasing is a Boolean which keeps track of whether
     the values in x_bound are increasing or decreasing.  This is used
     in the algorithm to decide exactly when the vertical line at x
     intersects the boundary; see below. */
  increasing=(x_bound[npts]>=x_bound[npts-1]);

  /* Circle the boundary checking for intersections. */
  for(i=0;i<npts;i++){

    /* An intersection is defined as x lying in the set
       (x_bound[i],x_bound[i+1]], or in the set
       [x_bound[i],x_bound[i+1]] if the direction of the curve is
       reversing itself about the value x_bound[i].  See the comments
       under test_cycle() for the conventions regarding direction
       reversal.  In this way, the return value of range() will always
       be even, and can be up to but no greater than the return value
       of test_cycle() applied to x_bound.

       Another way of looking at this: Vertical lines on the boundary
       are not treated as being truly vertical, but as sloping very
       slightly to the right as one proceeds towards increasing index
       values. */

    if((x>x_bound[i] && x<=x_bound[i+1]) ||
       (x<x_bound[i] && x>=x_bound[i+1]) ||
       (((x>=x_bound[i] && x<=x_bound[i+1]) ||
	 (x<=x_bound[i] && x>=x_bound[i+1])) &&
	increasing!=(x_bound[i+1]>=x_bound[i]))){

      /* When an intersection is found, interpolate to get the y
         value. */
      if(x_bound[i]==x_bound[i+1])
	y=y_bound[i];
      else
	y=y_bound[i]+(x-x_bound[i])*(y_bound[i+1]-y_bound[i])/
	  (x_bound[i+1]-x_bound[i]);
      if(n_int==0)
	*y_low=*y_high=y;
      else if(y<*y_low)
	*y_low=y;
      else if(y>*y_high)
	*y_high=y;
      n_int++;
    }
    increasing=(x_bound[i+1]>=x_bound[i]);
  }
  return n_int;
}


int get_height(int (*metric)(double , double , double *), double x,
	       double y, double dx, double *dy)

     /* This routine computes the height of a rectangular tile of
	specified width, such that the corners extend no more than
	unit proper distance from the centre, as measured by the
	metric function (*metric)().  get_height() returns 0 if
	successful, 1 if the requested location lies outside of the
	parameter space (leaving dy unchanged), or 2 if the requested
	width exceeds unit proper distance (in which case dy is set to
	0).

	The arguments are:

	metric(): Input.  This function computes the three independent
	  components of the distance metric matrix at a given point in
	  parameter space.  See the documentation for tiling_2d(),
	  above.

	x: Input.  The x-coordinate of the location in the
	  computational coordinate space.

	y: Input.  The y-coordinate of the location in the
	  computational coordinate space.

	dx: Input.  The requested width of the tile.

	dy: Output.  The height of the tile.

	The countour of unit proper distance is an ellipse defined by
	the metric equation:

	  1 = ds^2 = A*X*X + B*X*Y + C*Y*Y,

	where X,Y measure distances from the centre of the tile
	(whereas x,y measure the location of the tile); thus, dx=2*|X|
	and dy=2*|Y|.  A, B, and C are the metric coefficients
	computed by the metric() function.  In general for any
	(sufficiently small) X there will be two solutions for Y:

	  Y = (-B*X +/- sqrt(B*B*X*X-4*C*(A*X*X-1)))/(2*C).

	In order for the rectangle to be inscribed within the
	constant-match ellipse, these two solutions must both exist
	and have opposite sign; dy is then twice the magnitude of the
	smaller solution.  After a little work it can be deduced that:

	  dy = sqrt((B*dx/2C)^2 + (4-A*dx*dx)/C) - |B*dx/2C|,

	provided that this quantity is real and positive. */

{
  double m[3],temp,discriminant;

  /* Get the coefficients of the metric function.  Note that in the
     derivation above, A=m[0], B=m[1], C=m[2]. */
  if((*metric)(x,y,m))
    return 1;

  /* Compute the discriminant of the quadratic. */
  temp=fabs(0.5*m[1]*dx/m[2]);
  discriminant=(4.0-m[0]*dx*dx)/m[2]+temp*temp;
  if(discriminant<=0.0){
    *dy=0.0;
    return 2;
  }

  /* Compute the height of the tile. */
  if((*dy=sqrt(discriminant)-temp)<=0.0){
    *dy=0.0;
    return 2;
  }else
    return 0;
}


int get_width(int (*metric)(double , double , double *), double x,
	      double y, double *dx)

     /* This routine computes the width of the maximum-area
	rectangular tile lying entirely within unit proper distance of
	its centre, as defined by the distance function (*metric)().
	get_width() returns 0 if successful, or 1 if the requested
	location lies outside of the parameter space, or if the metric
	was otherwise undefined (in which case dx is left unchanged).

	The arguments are:

	metric(): Input.  This function computes the three independent
	  components of the distance metric matrix at a given point in
	  parameter space.  See the documentation for tiling_2d(),
	  above.

	x: Input.  The x-coordinate of the location in the
	  computational coordinate space.

	y: Input.  The y-coordinate of the location in the
	  computational coordinate space.

	dx: Output.  The width of the maximum-area rectangular tile.

	Essentially this routine finds the (positive) value of dx
	which maximizes the function dx*dy, where dy as a function of
	dx is given in the derivation of the get_height() routine,
	above:

	  dy = sqrt((B*dx/2C)^2 + (4-A*dx*dx)/C) - |B*dx/2C|.

	Requiring dx>0, we can rewrite this as:

	  dy = sqrt[ (|B/2C|^2 - A/C)*dx*dx + 4/C ] - |B/2C|*dx.

	From the fact that A and C are both positive, one can readily
	deduce that dy is positive for dx=0 and decreases to zero for
	some dx>0.  Thus dx*dy must have a positive local and absolute
	maximum for dx>0.  Setting the derivative of the area to zero,
	one gets after much algebra:

	  dx^2 = (2/A)/(1 +/- sqrt[B*B/4AC]).

	We note that B*B < 4AC (since the metric must be
	positive-definite), so each bracketed terms will always be
	positive.  It turns out that the larger root (- sign in the
	denominator) is spurious, so we have:

	  dx = sqrt[(2/A)/(1+sqrt[B*B/4AC])]. */

{
  double m[3],temp;

  /* Get the coefficients of the metric function.  Note that in the
     derivation above, A=m[0], B=m[1], C=m[2]. */
  if((*metric)(x,y,m))
    return 1;

  /* Compute the width of the tile. */
  temp=1.0+sqrt(m[1]*m[1]/(4.0*m[0]*m[2]));
  *dx=sqrt(2.0/(m[0]*temp));
  return 0;
}


int get_ellipse(int (*metric)(double , double , double *), double x,
		double y, double *r1, double *r2, double *theta)

     /* This routine computes the size and orientation of the
	elliptical contour at unit proper distance (as defined by the
	distance function metric()) about a specified location in
	parameter space.  This function returns 0 if successful, or 1
	if the requested location lies outside of the parameter space,
	or if the metric was undefined (leaving the output variables
	unchanged).

	The arguments are:

	metric(): Input.  This function computes the three independent
	  components of the distance metric matrix at a given point in
	  parameter space.  See the documentation for tiling_2d(),
	  above.

	x: Input.  The x-coordinate of the location in the
	  computational coordinate space.

	y: Input.  The y-coordinate of the location in the
	  computational coordinate space.

	r1: Output.  The semimajor axis of the ellipse.

	r2: Output.  The semiminor axis of the ellipse.

	theta: Output.  The angle counterclockwise from the x axis to
	  the semimajor axis. */

{
  double m[3],temp1,temp2,discriminant;

  /* Get the coefficients of the metric function. */
  if((*metric)(x,y,m))
    return 1;

  /* Compute the discriminant of the eigenvalues. */
  temp1=m[0]-m[2];
  discriminant=sqrt(temp1*temp1+m[1]*m[1]);

  /* Compute the semimajor and semiminor axis lengths. */
  temp2=m[0]+m[2];
  *r1=sqrt(2.0/(temp2-discriminant));
  *r2=sqrt(2.0/(temp2+discriminant));

  /* Compute the angle from the x-axis to the semimajor axis. */
  if(m[1]!=0.0)
    *theta=atan2(-m[1],discriminant-temp1);
  else if(m[0]<=m[2])
    *theta=0.0;
  else
    *theta=M_PI/2.0;
  return 0;
}


int get_plot_bounds(double *x_bound, double *y_bound, int npts,
		    double angle, double magnification,
		    double internal_margin, double *x_corner,
		    double *y_corner, int *horizontal_pages,
		    int *vertical_pages)

     /* This routine returns the number of pages of postscript output
        that will be generated by the routine plot_list(), below.

	The arguments are:

	x_bound: Input.  The array x_bound[0..npts] contains the x
	  components of a set of npts boundary points.  Note that the
	  array is of length npts+1; the [0] and [npts] index values
	  refer to the same point, so the array explicitly describes a
	  closed boundary; however, this is irrelevant to the current
	  routine.

	y_bound: Input.  The array y_bound[0..npts] contains the y
	  components of the boundary points, as above.

	npts: Input.  The number of points along the boundary.

	angle: Input.  The angle counterclockwise from the x-axis of
	  the parameter space to the horizontal axis of the plot.

	magnification: Input.  The scale factor of points (1/72 of an
	  inch) per unit coordinate distance in the parameter space.

	internal_margin: Input.  The extra space in the added around
	  the boundary, in the internal coordinate system.

	x_corner: Output.  The x-coordinate of the point which will
	  form the top left corner of the plot.

	y_corner: Output.  The y-coordinate of the point which will
	  form the top left corner of the plot.

	horizontal_pages: Output.  The number of pages in width of the
	  plot.

	vertical_pages: Output.  The number of pages in height of the
	  plot.

	Comments: Many of the important parameters of this function
	  are #defined constants at the top of the module.  The return
	  value of the function is equal to horizontal_pages times
	  vertical_pages. */

{
  double left,right,top,bottom,*x,*y,cosangle,sinangle,width,height,
    page_width,page_height;
  int i;

  /* Create local arrays of the rotated boundary curve. */
  x=(double *)malloc((npts+1)*sizeof(double));
  y=(double *)malloc((npts+1)*sizeof(double));
  cosangle=cos(angle);
  sinangle=sin(angle);
  for(i=0;i<=npts;i++){
    x[i]=x_bound[i]*cosangle-y_bound[i]*sinangle;
    y[i]=x_bound[i]*sinangle+y_bound[i]*cosangle;
  }

  /* Compute the number of pages required horizontally. */
  test_cycle(x,npts,&left,&right);
  width=magnification*(right-left+2.0*internal_margin);
  page_width=X_SIZE-2.0*X_MARGIN;
  *horizontal_pages=(int)(width/page_width)+1;

  /* Compute the number of pages required vertically. */
  test_cycle(y,npts,&bottom,&top);
  height=magnification*(top-bottom+2.0*internal_margin);
  page_height=Y_SIZE-2.0*Y_MARGIN;
  *vertical_pages=(int)(height/page_height)+1;

  /* Compute the x,y coordinates of the top left corner. */
  left-=internal_margin;
  top+=internal_margin;
  *x_corner=left*cosangle+top*sinangle;
  *y_corner=-left*sinangle+top*cosangle;

  /* De-allocate local arrays. */
  free(x);
  free(y);

  /* Return the total number of pages required. */
  return (*horizontal_pages)*(*vertical_pages);
}


int plot_list(double *x_bound, double *y_bound, int npts,
	      struct tile *head, int n_tiles, double angle,
	      double magnification, int plot_boundary,
	      int plot_tiles, int plot_ellipses, int plot_flags,
	      const char *psfile)

     /* This routine generates a postscript file displaying a
	parameter space, the brickwork mesh of tiles covering it, and
	an overlapping mesh of elliptical contours of unit proper
	radius, circumscribing each tile.  The function returns the
	number of pages of postscript output.

	The arguments are:

	x_bound: Input.  The array x_bound[0..npts] contains the x
	  components of a set of npts boundary points.  Note that the
	  array is of length npts+1; the [0] and [npts] index values
	  refer to the same point, so the array explicitly describes a
	  closed boundary; however, this is irrelevant to the current
	  routine.

	y_bound: Input.  The array y_bound[0..npts] contains the y
	  components of the boundary points, as above.

	npts: Input.  The number of points along the boundary.

	head: Input.  The head of the linked list of tiles to be
	  plotted.

	n_tiles: Input.  The number of tiles to be plotted from the
	  list.

	angle: Input.  The angle counterclockwise from the x-axis of
	  the parameter space to the horizontal axis of the plot.

	magnification: Input.  The scale factor of points (1/72 of an
	  inch) per unit coordinate distance in the parameter space.

	plot_boundary: Input.  1 if boundary is to be shown, 0
	  otherwise.

	plot_tiles: Input.  1 if rectangular brickwork of tiles is to
	  be shown, 0 otherwise.

	plot_ellipses: Input.  1 if the overlapping match contours are
	  to be shown, 0 otherwise.

	plot_flags: Input.  If nonzero, indicates the size of dot used
	  to mark flagged tiles (in points = 1/72 inches).  If zero,
	  flags are ignored.

	psfile: Input.  The name of the postscript file created.

	Comments: Many of the important parameters of this function
	  are #defined constants at the top of the module. */

{
  double internal_margin,x_corner,y_corner;
  int pages,horizontal_pages,vertical_pages,n,i,j,k;
  int nobj_boundary=0,nobj_tiles=0,nobj_ellipses=0,nobj_flags=0;
  time_t tp;
  struct tile *here;
  FILE *fpout;

  /* Deternime the positioning parameters of the plot. */
  internal_margin=0.0;
  if(plot_tiles || plot_ellipses)
    for(here=head,i=0;here!=NULL && i<n_tiles;here=(*here).next,i++)
      if((*here).r1>internal_margin)
	internal_margin=(*here).r1;
  pages=get_plot_bounds(x_bound,y_bound,npts,angle,magnification,
			internal_margin,&x_corner,&y_corner,
			&horizontal_pages,&vertical_pages);

  /* Open postscript file and write postscript header. */
  fpout=fopen(psfile,"w");
  fprintf(fpout,"%%!PS-Adobe-2.0\n");
  fprintf(fpout,"%%%%Creator: plot_list()\n");
  fprintf(fpout,"%%%%Title: ");
  fprintf(fpout,psfile);
  fprintf(fpout,"\n");
  fprintf(fpout,"%%%%Creation Date: ");
  time(&tp);
  fprintf(fpout,ctime(&tp));
  fprintf(fpout,"\n");
  fprintf(fpout,"%%%%Pages: %i (%i by %i)\n",pages,horizontal_pages,
	  vertical_pages);
  fprintf(fpout,"%%%%Pageorder: Ascent\n");
  fprintf(fpout,"%%%%BoundingBox: 0 0 %i %i\n",X_SIZE,Y_SIZE);
  fprintf(fpout,"%%%%EndComments\n\n");

  /* Set line width to zero. */
  fprintf(fpout,"0 setlinewidth\n");

  /* Define and save a clipping box around the plot region. */
  fprintf(fpout,"newpath\n");
  fprintf(fpout,"%i %i moveto\n",X_MARGIN,Y_MARGIN);
  fprintf(fpout,"%i %i lineto\n",X_SIZE-X_MARGIN,Y_MARGIN);
  fprintf(fpout,"%i %i lineto\n",X_SIZE-X_MARGIN,Y_SIZE-Y_MARGIN);
  fprintf(fpout,"%i %i lineto\n",X_MARGIN,Y_SIZE-Y_MARGIN);
  fprintf(fpout,"closepath clip gsave\n\n");

  if(plot_boundary){
    /* Define macros to plot the boundary of the space. */
    nobj_boundary=0;
    fprintf(fpout,"/boundary%i {newpath\n",nobj_boundary);
    fprintf(fpout,"%f %f moveto\n",x_bound[0],y_bound[0]);
    for(i=1,j=8;i<=npts;i++,j+=3){
      if(j>N_OBJ){
	fprintf(fpout,"stroke} def\n\n");
	nobj_boundary++;
	fprintf(fpout,"/boundary%i {newpath\n",nobj_boundary);
	fprintf(fpout,"%f %f moveto\n",x_bound[i-1],y_bound[i-1]);
	j=8;
      }
      fprintf(fpout,"%f %f lineto\n",x_bound[i],y_bound[i]);
    }
    fprintf(fpout,"stroke} def\n\n");
    nobj_boundary++;
  }

  if(plot_tiles){
    /* Define a macro to plot a single rectangular tile.
       Usage: (x1) (y1) (x2) (y2) box
       where: (x1),(y1) = coordinates of one corner of the box,
              (x1),(y1) = coordinates of the opposite corner. */
    fprintf(fpout,"/box {newpath 1 index 3 index moveto 2 copy lineto"
	    " 3 index exch lineto pop lineto closepath stroke}"
	    " def\n\n");

    /* Define macros to plot the mesh of tiles. */
    nobj_tiles=0;
    fprintf(fpout,"/tiles%i {\n",nobj_tiles);
    for(here=head,i=0,j=5;here!=NULL && i<n_tiles;
	here=(*here).next,i++,j+=5){
      if(j>N_OBJ){
	fprintf(fpout,"} def\n\n");
	nobj_tiles++;
	fprintf(fpout,"/tiles%i {\n",nobj_tiles);
	j=5;
      }
      fprintf(fpout,"%f %f %f %f box\n",
	      (*here).x-0.5*(*here).dx,(*here).y-0.5*(*here).dy,
	      (*here).x+0.5*(*here).dx,(*here).y+0.5*(*here).dy);
    }
    fprintf(fpout,"} def\n\n");
    nobj_tiles++;
  }

  if(plot_ellipses){
    /* Define a macro to plot a single elliptical contour.
       Usage: (r1) (r2) (theta) (x) (y) el
       where: (r1),(r2) = semimajor and semiminor axis lengths,
              (theta)   = angle from x to semimajor axis,
              (x),(y)   = coordinates of centre. */
    fprintf(fpout,"/el {gsave translate rotate scale newpath 0 0 1 0"
	    " 360 arc stroke grestore} def\n\n");

    /* Define macros to plot the mesh of overlapping contours. */
    nobj_ellipses=0;
    fprintf(fpout,"/ellipses%i {\n",nobj_ellipses);
    for(here=head,i=0,j=6;here!=NULL && i<n_tiles;
	here=(*here).next,i++,j+=6){
      if(j>N_OBJ){
	fprintf(fpout,"} def\n\n");
	nobj_ellipses++;
	fprintf(fpout,"/ellipses%i {\n",nobj_ellipses);
	j=6;
      }
      fprintf(fpout,"%f %f %f %f %f el\n",(*here).r1,(*here).r2,
	      RtD*(*here).theta,(*here).x,(*here).y);
    }
    fprintf(fpout,"} def\n\n");
    nobj_ellipses++;
  }

  if(plot_flags){
    /* Define a macro to plot a single flag indicator.
       Usage: (x) (y) dot
       where: (x),(y) = coordinates of the flagged tile. */
    fprintf(fpout,"/dot {newpath %f 0 360 arc closepath fill}"
	    " def\n\n",0.5*plot_flags/magnification);

    /* Define macros to plot indicators on all flagged tiles. */
    nobj_flags=0;
    fprintf(fpout,"/flags%i {\n",nobj_flags);
    for(here=head,i=0,j=3;here!=NULL && i<n_tiles;
	here=(*here).next,i++)
      if((*here).flag){
	if(j>N_OBJ){
	  fprintf(fpout,"} def\n\n");
	  nobj_flags++;
	  fprintf(fpout,"/flags%i {\n",nobj_flags);
	  j=3;
	}
	fprintf(fpout,"%f %f dot\n",(*here).x,(*here).y);
	j+=3;
      }
    fprintf(fpout,"} def\n\n");
    nobj_flags++;
  }

  /* Print each page: */
  for(i=0,n=1;i<vertical_pages;i++)
    for(j=0;j<horizontal_pages;j++,n++){
      fprintf(fpout,"%%%%Page: graph %i\n\n",n);

      /* Locate origin (top left corner) of the plot space relative to
         the current page, rescale, and rotate. */
      fprintf(fpout,"grestore gsave\n");
      fprintf(fpout,"%f %f translate\n",
	      X_MARGIN-j*(X_SIZE-2.0*X_MARGIN),
	      Y_SIZE-Y_MARGIN+i*(Y_SIZE-2.0*Y_MARGIN));
      fprintf(fpout,"%f %f scale\n",magnification,magnification);
      fprintf(fpout,"%f rotate\n",RtD*angle);
      fprintf(fpout,"%f %f translate\n",-x_corner,-y_corner);

      /* Plot the boundary. */
      if(plot_boundary)
	for(k=0;k<nobj_boundary;k++)
	  fprintf(fpout,"boundary%i\n",k);

      /* Plot the mesh of tiles. */
      if(plot_tiles)
	for(k=0;k<nobj_tiles;k++)
	  fprintf(fpout,"tiles%i\n",k);

      /* Plot the overlapping ellipses. */
      if(plot_ellipses)
	for(k=0;k<nobj_ellipses;k++)
	  fprintf(fpout,"ellipses%i\n",k);

      /* Plot the flag indicators. */
      if(plot_flags)
	for(k=0;k<nobj_flags;k++)
	  fprintf(fpout,"flags%i\n",k);

      /* Finish the page. */
      fprintf(fpout,"showpage\n\n");
    }

  /* Close the file and return the number of pages printed. */
  fclose(fpout);
  return --n;
}
