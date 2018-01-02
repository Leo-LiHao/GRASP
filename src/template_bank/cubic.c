/* GRASP: Copyright 1997,1998  Bruce Allen */
static char *rcsid="$Id: cubic.c,v 1.8 1998/02/04 21:18:18 ballen Exp $\
\n$Name: RELEASE_1_9_8 $";

#include "grasp.h"

void generate_cubic(struct cubic_grid grid, char *detectors_file,
		    const char *outfile, const char *logfile)

     /* This routine computes the coefficients of the cubic fit to the
        match function on a mesh in parameter space, and writes the
        results to an ASCII textfile, suitable for reading by the
        routine read_cubic() (below).

	The arguments are:

	grid: Input/Output.  This structure contains the parameters
	  for the computation of the cubic fits.  All of the fields
	  except for grid.dm and grid.coef must be set; those fields
	  are the ones that are computed.

	detectors_file: Input.  The name of a data file containing
	  detector site information, such as detectors.dat.  This is
	  used to get a noise file for computing the match function.

	outfile: Input.  The name of the output file to which the
	  coefficients and related information will be written.

	logfile: Input.  The name of a log file which tracks the
	  progress of this routine (since it can take several hours to
	  generate a reasonable grid).

	The output file is an ASCII textfile containing the fields of
	the structure grid.  Each field except the coef field is
	printed on a separate line of the output file.  The coef data
	is written as 0.5*(grid.n)*(grid.n+1) lines of 10 floating
	point numbers; each line represents the 10 coefficients
	coef[i][j][0..9] for a given i,j.  The lines are ordered by
	increasing j from 0..i for each i from 0..grid.n-1.  Integers
	are printed exactly; floats are printed in 10-digit precision
	exponential notation.

	One should also note that this routine can take quite a long
	time to run: on a 100 MHz pentium it typically takes 10 to 15
	minutes per point in the grid, and there are
	(grid.n)*(grid.n+1)/2 grid points in the triangular parameter
	space.  This is the reason for creating the log file to track
	the routine's progress. */

{
  float m1,m2,site_parameters[9],c,s,co[10];
  int i,j,k;
  char noise_file[128],whiten_file[128],site_name[128];
  FILE *fplog,*fpout;

  /* Open output and log files. */
  fpout=fopen(outfile,"w");
  fplog=fopen(logfile,"w");

  /* Compute the mass spacing for the grid. */
  grid.dm=(grid.m_mx-grid.m_mn)/(grid.n-1);

  /* This completes the grid structure except for the actual
     coefficients, so write the grid structure to the output file. */
  fprintf(fpout,"struct cubic_grid data file, created by"
	  " generate_cubic().\n");
  fprintf(fpout,"(Manual editing may cause corruption.)\n");
  fprintf(fpout,"n: %i\n",grid.n);
  fprintf(fpout,"m_mn: %15.9e\n",grid.m_mn);
  fprintf(fpout,"m_mx: %15.9e\n",grid.m_mx);
  fprintf(fpout,"dm: %15.9e\n",grid.dm);
  fprintf(fpout,"match: %15.9e\n",grid.match);
  fprintf(fpout,"angle: %15.9e\n",grid.angle);
  fprintf(fpout,"order: %i\n",grid.order);
  fprintf(fpout,"srate: %15.9e\n",grid.srate);
  fprintf(fpout,"flo: %15.9e\n",grid.flo);
  fprintf(fpout,"ftau: %15.9e\n",grid.ftau);
  fprintf(fpout,"detector: %i\n",grid.detector);
  fprintf(fpout,"coef[i=0..%i][j=0..i][k=0..9]:\n",grid.n-1);
  fflush(fpout);

  /* Now get the name of the noise file for computing the match
     function. */
  detector_site(detectors_file,grid.detector,site_parameters,
		site_name,noise_file,whiten_file);
  fprintf(fplog,"Using %s noise curve from the file %s.\n",site_name,
	  noise_file);
  fflush(fplog);

  /* Allocate memory for the coefficient grid. */
  grid.coef=(float ***)malloc(grid.n*sizeof(float **));
  for(i=0;i<grid.n;i++){
    grid.coef[i]=(float **)malloc(grid.n*sizeof(float *));
    for(j=0;j<=i;j++)
      grid.coef[i][j]=(float *)malloc(10*sizeof(float));
  }

  /* Compute the rotation transformation. */
  c=cos(grid.angle);
  s=sin(grid.angle);

  /* Now compute the coefficients, and write them to the file. */
  fprintf(fplog,"Generating match coefficients at %i points:\n",
	  grid.n*(grid.n+1)/2);
  fflush(fplog);
  for(m1=grid.m_mn,i=0;i<grid.n;m1+=grid.dm,i++){
    for(m2=grid.m_mn,j=0;j<=i;m2+=grid.dm,j++){

      /* Find the unrotated coefficients. */
      match_cubic(m1,m2,grid.match,grid.order,grid.srate,grid.flo,
		  grid.ftau,noise_file,co+7,co+8,co+9,co);

      /* Rotate quadratic components. */
      grid.coef[i][j][0]=co[0]*c*c+co[1]*c*s+co[2]*s*s;
      grid.coef[i][j][1]=2.0*(co[2]-co[0])*c*s+co[1]*(c*c-s*s);
      grid.coef[i][j][2]=co[0]*s*s-co[1]*c*s+co[2]*c*c;

      /* Rotate cubic components. */
      grid.coef[i][j][3]=co[3]*c*c*c+co[4]*s*s*s+co[5]*c*c*s
	+co[6]*c*s*s;
      grid.coef[i][j][4]=-co[3]*s*s*s+co[4]*c*c*c+co[5]*c*s*s
	-co[6]*c*c*s;
      grid.coef[i][j][5]=(2.0*co[6]-3.0*co[3])*c*c*s
	+(3.0*co[4]-2.0*co[5])*c*s*s+co[5]*c*c*c-co[6]*s*s*s;
      grid.coef[i][j][6]=(3.0*co[3]-2.0*co[6])*c*s*s
	+(3.0*co[4]-2.0*co[5])*c*c*s+co[5]*s*s*s+co[6]*c*c*c;

      /* Compute ellipse. */
      grid.coef[i][j][7]=co[7];
      grid.coef[i][j][8]=co[8];
      grid.coef[i][j][9]=co[9]-grid.angle;

      /* Write coefficients to data file. */
      fprintf(fpout,"[%i][%i]: ",i,j);
      for(k=0;k<10;k++)
	fprintf(fpout,"%15.9e ",grid.coef[i][j][k]);
      fprintf(fpout,"\n");
      fflush(fpout);
      fprintf(fplog,".");
      fflush(fplog);
    }
    fprintf(fplog,"\n");
    fflush(fplog);
  }
  fprintf(fplog,"\n");
  fflush(fplog);

  /* Free the coefficient array, close files, and exit. */
  free_cubic(grid);
  fclose(fpout);
  fclose(fplog);
  return;
}


int regenerate_cubic(char *detectors_file, const char *infile,
		     const char *outfile, const char *logfile)

     /* It may happen that the routine generate_cubic() terminates
	before completing the entire coefficient grid.  Since each
	grid point takes so long to compute, it would be foolish to
	discard those already generated.  This routine, therefore,
	reads in a partially-complete data file, and then continues
	the computation where generate_cubic() left off.  The results
	are written to a new data file (leaving the original file
	incomplete).  It returns 0 upon successful completion, or 1 if
	the data file was absent or corrupt.

	regenerate_cubic() also creates its own log file to track its
	progress.

	The arguments are:

	detectors_file: Input.  The name of a data file containing
	  detector site information, such as detectors.dat.  This is
	  used to get a noise file for computing the match function.

	infile: Input.  The name of the incomplete data file of
	  coefficients.

	outfile: Input.  The name of the data file where this routine
	  stores its results.

	logfile: Input.  The name of a log file which tracks the
	  progress of this routine. */

{
  struct cubic_grid grid;
  float m1,m2,site_parameters[9],c,s,co[10];
  int i,j,k,n,n_start,flag;
  char noise_file[128],whiten_file[128],site_name[128],checkstr[128],
    readstr[128];
  FILE *fpin,*fpout,*fplog;

  /* Open files. */
  if((fpin=fopen(infile,"r"))==NULL)
    return 1;
  fpout=fopen(outfile,"w");
  fplog=fopen(logfile,"w");

  /* Read the coefficient grid parameters from the input file, and
     write them to the output file. */
  if(fscanf(fpin,
	    "struct cubic_grid data file, created by %*s\n"
	    "(Manual editing may cause corruption.)\n"
	    "n: %i\n"
	    "m_mn: %f\n"
	    "m_mx: %f\n"
	    "dm: %f\n"
	    "match: %f\n"
	    "angle: %f\n"
	    "order: %i\n"
	    "srate: %f\n"
	    "flo: %f\n"
	    "ftau: %f\n"
	    "detector: %i\n",&grid.n,&grid.m_mn,&grid.m_mx,&grid.dm,
	    &grid.match,&grid.angle,&grid.order,&grid.srate,&grid.flo,
	    &grid.ftau,&grid.detector)!=11)
    return 1;
  sprintf(checkstr,"coef[i=0..%i][j=0..i][k=0..9]:\n",grid.n-1);
  fscanf(fpin,"%127s",readstr);
  if(strcmp(checkstr,readstr))
    return 1;

  fprintf(fpout,"struct cubic_grid data file, created by"
	  " regenerate_cubic().\n");
  fprintf(fpout,"(Manual editing may cause corruption.)\n");
  fprintf(fpout,"n: %i\n",grid.n);
  fprintf(fpout,"m_mn: %15.9e\n",grid.m_mn);
  fprintf(fpout,"m_mx: %15.9e\n",grid.m_mx);
  fprintf(fpout,"dm: %15.9e\n",grid.dm);
  fprintf(fpout,"match: %15.9e\n",grid.match);
  fprintf(fpout,"angle: %15.9e\n",grid.angle);
  fprintf(fpout,"order: %i\n",grid.order);
  fprintf(fpout,"srate: %15.9e\n",grid.srate);
  fprintf(fpout,"flo: %15.9e\n",grid.flo);
  fprintf(fpout,"ftau: %15.9e\n",grid.ftau);
  fprintf(fpout,"detector: %i\n",grid.detector);
  fprintf(fpout,"coef[i=0..%i][j=0..i][k=0..9]:\n",grid.n-1);
  fflush(fpout);

  /* Get the name of the noise file for computing the match
     function. */
  detector_site(detectors_file,grid.detector,site_parameters,
		site_name,noise_file,whiten_file);
  fprintf(fplog,"Using %s noise curve from the file %s.\n",site_name,
	  noise_file);
  fflush(fplog);

  /* Allocate memory for the coefficient grid. */
  grid.coef=(float ***)malloc(grid.n*sizeof(float **));
  for(i=0;i<grid.n;i++){
    grid.coef[i]=(float **)malloc(grid.n*sizeof(float *));
    for(j=0;j<=i;j++)
      grid.coef[i][j]=(float *)malloc(10*sizeof(float));
  }

  /* Compute the rotation transformation. */
  c=cos(grid.angle);
  s=sin(grid.angle);

  /* Read the existing coefficients from the input file. */
  n_start=0.5*grid.n*(grid.n+1)+1;
  for(i=0,n=0,flag=1;i<grid.n && flag;i++)
    for(j=0;j<=i && flag;j++,n++){
      sprintf(checkstr,"[%i][%i]:",i,j);
      fscanf(fpin,"%127s",readstr);
      if(strcmp(checkstr,readstr)){
	n_start=n;
	flag=0;
      }
      for(k=0;k<10;k++)
	if(fscanf(fpin,"%f",grid.coef[i][j]+k)!=1){
	  n_start=n;
	  flag=0;
	}
    }

  /* Write the existing coefficients to the output file, and compute
     and write the remaining coefficients. */
  fprintf(fplog,"Generating match coefficients at %i points:\n",
	  grid.n*(grid.n+1)/2);
  fflush(fplog);
  for(m1=grid.m_mn,i=0,n=0;i<grid.n;m1+=grid.dm,i++){
    for(m2=grid.m_mn,j=0;j<=i;m2+=grid.dm,j++,n++){

      /* Don't recompute components which have been read in. */
      if(n>=n_start){

	/* Find the unrotated coefficients. */
	match_cubic(m1,m2,grid.match,grid.order,grid.srate,grid.flo,
		    grid.ftau,noise_file,co+7,co+8,co+9,co);

	/* Rotate quadratic components. */
	grid.coef[i][j][0]=co[0]*c*c+co[1]*c*s+co[2]*s*s;
	grid.coef[i][j][1]=2.0*(co[2]-co[0])*c*s+co[1]*(c*c-s*s);
	grid.coef[i][j][2]=co[0]*s*s-co[1]*c*s+co[2]*c*c;

	/* Rotate cubic components. */
	grid.coef[i][j][3]=co[3]*c*c*c+co[4]*s*s*s+co[5]*c*c*s
	  +co[6]*c*s*s;
	grid.coef[i][j][4]=-co[3]*s*s*s+co[4]*c*c*c+co[5]*c*s*s
	  -co[6]*c*c*s;
	grid.coef[i][j][5]=(2.0*co[6]-3.0*co[3])*c*c*s
	  +(3.0*co[4]-2.0*co[5])*c*s*s+co[5]*c*c*c-co[6]*s*s*s;
	grid.coef[i][j][6]=(3.0*co[3]-2.0*co[6])*c*s*s
	  +(3.0*co[4]-2.0*co[5])*c*c*s+co[5]*s*s*s+co[6]*c*c*c;

	/* Compute ellipse. */
	grid.coef[i][j][7]=co[7];
	grid.coef[i][j][8]=co[8];
	grid.coef[i][j][9]=co[9]-grid.angle;
      }

      /* Write coefficients to data file. */
      fprintf(fpout,"[%i][%i]: ",i,j);
      for(k=0;k<10;k++)
	fprintf(fpout,"%15.9e ",grid.coef[i][j][k]);
      fprintf(fpout,"\n");
      fflush(fpout);
      fprintf(fplog,".");
      fflush(fplog);
    }
    fprintf(fplog,"\n");
    fflush(fplog);
  }
  fprintf(fplog,"\n");
  fflush(fplog);

  /* Free the coefficient array, close files, and exit. */
  free_cubic(grid);
  fclose(fpin);
  fclose(fpout);
  fclose(fplog);
  return 0;
}


int read_cubic(struct cubic_grid *grid, const char *infile)

     /* This routine reads a textfile generated by generate_cubic() or
	write_cubic(), stores the data in a variable grid of type
	struct cubic_grid, and passes this structure back.
	read_cubic() itself returns 0 after successful completion, or
	1 if the data file was absent or corrupt (in which case grid
	is left unchanged).

	Note that memory for the coefficient array grid.coef is
	allocated in this routine; to free this memory, call
	free_cubic(grid).  Allocating and de-allocating this array
	requires some care: since grid.coef[i][j][k] is necessarily
	symmetric in the first two indecies, the pointers
	grid.coef[i][j] and grid.coef[j][i] have been explicitly set
	to point to the same memory location, in order to save memory.

	The arguments are:

	grid: Output.  The coefficient array and related data read
	  from the data file.

	infile: Input.  The name of the data file. */

{
  struct cubic_grid temp;
  int i,j,k;
  char checkstr[128],readstr[128];
  FILE *fpin;

  /* Open data file and read the coefficient grid parameters. */
  if((fpin=fopen(infile,"r"))==NULL)
    return 1;
  if(fscanf(fpin,
	    "struct cubic_grid data file, created by %*s\n"
	    "(Manual editing may cause corruption.)\n"
	    "n: %i\n"
	    "m_mn: %f\n"
	    "m_mx: %f\n"
	    "dm: %f\n"
	    "match: %f\n"
	    "angle: %f\n"
	    "order: %i\n"
	    "srate: %f\n"
	    "flo: %f\n"
	    "ftau: %f\n"
	    "detector: %i\n",&temp.n,&temp.m_mn,&temp.m_mx,&temp.dm,
	    &temp.match,&temp.angle,&temp.order,&temp.srate,&temp.flo,
	    &temp.ftau,&temp.detector)!=11)
    return 1;
  sprintf(checkstr,"coef[i=0..%i][j=0..i][k=0..9]:",temp.n-1);
  fscanf(fpin,"%127s",readstr);
  if(strcmp(checkstr,readstr))
    return 1;

  /* Allocate memory for the coefficient grid. */
  temp.coef=(float ***)malloc(temp.n*sizeof(float **));
  for(i=0;i<temp.n;i++){
    temp.coef[i]=(float **)malloc(temp.n*sizeof(float *));
    for(j=0;j<=i;j++)
      temp.coef[i][j]=temp.coef[j][i]
	=(float *)malloc(10*sizeof(float));
  }

  /* Read the coefficients from the data file. */
  for(i=0;i<temp.n;i++)
    for(j=0;j<=i;j++){
      sprintf(checkstr,"[%i][%i]:",i,j);
      fscanf(fpin,"%127s",readstr);
      if(strcmp(checkstr,readstr)){
	free_cubic(temp);
	return 1;
      }
      for(k=0;k<10;k++){
	if(fscanf(fpin,"%f",temp.coef[i][j]+k)!=1){
	  free_cubic(temp);
	  return 1;
	}
      }
    }

  /* Close the file and return the grid structure. */
  fclose(fpin);
  *grid=temp;
  return 0;
}


void write_cubic(struct cubic_grid grid, const char *outfile)

     /* This routine takes an existing cubic coefficient grid stored
	in a structure of type struct cubic_grid, and writes it to an
	ASCII textfile in a format suitable for reading by
	read_cubic().

	The arguments are:

	grid: Input.  The coefficient array and related data to be
	  written to the data file.

	outfile: Input.  The name of the data file. */

{
  int i,j,k;
  FILE *fpout;

  /* Open data file and write the coefficient grid parameters. */
  fpout=fopen(outfile,"w");
  fprintf(fpout,"struct cubic_grid data file, created by"
	  " write_cubic().\n");
  fprintf(fpout,"(Manual editing may cause corruption.)\n");
  fprintf(fpout,"n: %i\n",grid.n);
  fprintf(fpout,"m_mn: %15.9e\n",grid.m_mn);
  fprintf(fpout,"m_mx: %15.9e\n",grid.m_mx);
  fprintf(fpout,"dm: %15.9e\n",grid.dm);
  fprintf(fpout,"match: %15.9e\n",grid.match);
  fprintf(fpout,"angle: %15.9e\n",grid.angle);
  fprintf(fpout,"order: %i\n",grid.order);
  fprintf(fpout,"srate: %15.9e\n",grid.srate);
  fprintf(fpout,"flo: %15.9e\n",grid.flo);
  fprintf(fpout,"ftau: %15.9e\n",grid.ftau);
  fprintf(fpout,"detector: %i\n",grid.detector);
  fprintf(fpout,"coef[i=0..%i][j=0..i][k=0..9]:\n",grid.n-1);

  /* Write the coefficients themselves. */
  for(i=0;i<grid.n;i++)
    for(j=0;j<=i;j++){
      fprintf(fpout,"[%i][%i]: ",i,j);
      for(k=0;k<10;k++)
	fprintf(fpout,"%15.9e ",grid.coef[i][j][k]);
      fprintf(fpout,"\n");
    }

  /* Close the file. */
  fclose(fpout);
  return;
}


int get_cubic(float m1, float m2, struct cubic_grid grid, float *coef)

     /* This routine computes the coefficients of a cubic fit to the
        match function at a specified point in parameter space, by
        linear interpolation of precomputed coefficients on a grid in
        parameter space.  It returns 0 if successfully executed, or 1
        if the point (m1,m2) lies outside of the grid.  In the latter
        case, get_cubic() will compute extrapolated coefficients, but
        these are unreliable.

	The arguments are:

	m1: Input.  One of the binary mass coordinates of the
	  requested point in parameter space.

	m2: Input.  The other binary mass coordinates of the requested
	  point in parameter space.

	grid: Input.  The data structure containing the precomputed
	  coefficients and the information required to retrieve them.

	coef: Output.  The array coef[0..9] is filled with the
	  interpolated coefficients. */

{
  float x,y,temp1,temp2;
  int i,j,k;

  /* Get the location of the requested point in grid coordinates. */
  x=(m1-grid.m_mn)/grid.dm;
  y=(m2-grid.m_mn)/grid.dm;

  /* Get the indecies of the lower-left corner of the grid patch
     containing the requested point (or closest to it if it lies
     outside the grid). */
  if(x<0)
    i=0;
  else if(x>=grid.n-1)
    i=grid.n-2;
  else
    i=(int)x;
  if(y<0)
    j=0;
  else if(y>=grid.n-1)
    j=grid.n-2;
  else
    j=(int)y;

  /* Get the location of the requested point with respect to the lower
     left corner of this patch. */
  x-=i;
  y-=j;

  /* Interpolate (or extrapolate) the four corners of the patch to get
     the coefficients at the requested point. */
  for(k=0;k<10;k++){
    temp1=(1.0-x)*grid.coef[i][j][k]+x*grid.coef[i+1][j][k];
    temp2=(1.0-x)*grid.coef[i][j+1][k]+x*grid.coef[i+1][j+1][k];
    coef[k]=(1.0-y)*temp1+y*temp2;
  }

  /* Return whether or not the requested point was inside the grid. */
  return (x<0.0 || x>1.0 || y<0.0 || y>1.0);
}


void free_cubic(struct cubic_grid grid)

     /* Frees the memory allocated to the array grid.coef by
        read_cubic().

	The argument is:

	grid: Input.  The cubic_grid structure whose coefficient array
	  is to be freed. */

{
  int i,j;

  for(i=0;i<grid.n;i++){
    for(j=0;j<=i;j++)
      free(grid.coef[i][j]);
    free(grid.coef[i]);
  }
  free(grid.coef);
  return;
}


void transform_cubic(struct cubic_grid *grid, float angle,
		     float match)

     /* This routine applies a rotation to the coefficients stored in
	(*grid), and rescales the equimatch ellipses to a new match
	level.

	The arguments are:

	grid: Input/Output.  The structure containing the coefficients
	  to be transformed.

	angle: Input.  The new value of grid.angle (the elements
	  (*grid).coef[i][j][0..6,9] will be transformed to fit this
	  new angle).

	match: Input.  The new value of grid.match (the elements
	  (*grid).coef[i][j][7,8] will be rescaled according to this
	  value).

	Comments: The result of this transformation is not quite the
	  same as if the grid were originally generated with the new
	  value of the match.  During initial generation of the grid,
	  the match field also specifies the domain over which the
	  cubic fit is made, as well as setting the scale for the
	  equimatch ellipse axes.  This routine only rescales the
	  ellipses; it does not regenerate a cubic fit over a new
	  match range. */

{
  float c,s,scale,*co,coef[10];
  int i,j,k;

  /* Compute the relative rotation angle. */
  angle-=(*grid).angle;
  (*grid).angle+=angle;

  /* Compute the components of the rotation transformation.  These
     temporary variables are used primarily to simplify the notation;
     I leave it to the compiler to optimize the computations. */
  c=cos(angle);
  s=sin(angle);

  /* Compute the ellipse scaling factor. */
  scale=pow((1.0-match)/(1.0-(*grid).match),0.5);
  (*grid).match=match;

  /* Now apply the transformation over the entire grid: */
  for(i=0;i<(*grid).n;i++)
    for(j=0;j<=i;j++){
      co=(*grid).coef[i][j];

      /* Transform quadratic components. */
      coef[0]=co[0]*c*c+co[1]*c*s+co[2]*s*s;
      coef[1]=2.0*(co[2]-co[0])*c*s+co[1]*(c*c-s*s);
      coef[2]=co[0]*s*s-co[1]*c*s+co[2]*c*c;

      /* Transform cubic components. */
      coef[3]=co[3]*c*c*c+co[4]*s*s*s+co[5]*c*c*s+co[6]*c*s*s;
      coef[4]=-co[3]*s*s*s+co[4]*c*c*c+co[5]*c*s*s-co[6]*c*c*s;
      coef[5]=(2.0*co[6]-3.0*co[3])*c*c*s+(3.0*co[4]-2.0*co[5])*c*s*s
	+co[5]*c*c*c-co[6]*s*s*s;
      coef[6]=(3.0*co[3]-2.0*co[6])*c*s*s+(3.0*co[4]-2.0*co[5])*c*c*s
	+co[5]*s*s*s+co[6]*c*c*c;

      /* Rescale and rotate ellipses. */
      coef[7]=co[7]*scale;
      coef[8]=co[8]*scale;
      coef[9]=co[9]-angle;

      /* Store new coefficient values. */
      for(k=0;k<10;k++)
	co[k]=coef[k];
    }
  return;
}
