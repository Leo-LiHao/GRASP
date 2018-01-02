/* GRASP: Copyright 1997,1998  Bruce Allen */
static char *rcsid="$Id: postscript.c,v 1.2 1998/01/23 17:59:51 ballen Exp $\
\n$Name: RELEASE_1_9_8 $";

#include "grasp.h"

#define PS_LM 36. /* left margin, units 1/72 inch */
#define PS_RM 36. /* right margin, units 1/72 inch */
#define PS_BM 72. /* bottom margin, units 1/72 inch */
#define PS_PW 540. /* printable page width, units 1/72 inch */
#define PS_CH 10000. /* clip height, units 1/72 inch */
#define PS_OX 100. /* x origin on page, units 1/72 inch */
#define PS_OY 396. /* y origin on page, units 1/72 inch */



void plot_template(char *filename,struct Scope Grid,int npages,int number) {
int i;
FILE *fps;
char *date;
time_t timer;
float magnification;
void plot_template_postscript(FILE *fp,struct Scope,int,float,int,int,char*,char*);

time(&timer);
date=ctime(&timer);
fps=fopen(filename,"w");

/* create the output file, printing a warning if there is a problem, and giving up */
if (fps==NULL) {
   GR_start_error("plot_template()",rcsid,__FILE__,__LINE__);
   GR_report_error("Unable to open file %s for postscript output.\n",filename);
   GR_report_error("Giving up on outputing postscript, but continuing execution...\n");
   GR_end_error();
   return;
}
fprintf(fps,"%%!PS-Adobe-2.0\n");
fprintf(fps,"%%%%Creator: plot_template()\n");
fprintf(fps,"%%%%Title: %s\n",filename);
fprintf(fps,"%%%%Creation Date: %s\n",date);
fprintf(fps,"%%%%Pages: %d\n",npages);
fprintf(fps,"%%%%Pageorder: Ascent\n");
fprintf(fps,"%%%%BoundingBox: 0 0 612 792\n");
fprintf(fps,"%%%%EndComments\n");


/* set magnification */
if (npages>1)
	magnification=npages-1.0;
else
	magnification=1.0;

/* create the npages of output */
for (i=1;i<=npages;i++)
	plot_template_postscript(fps,Grid,1,magnification,i,number,filename,date);

/* close the output file, and return */
GR_start_error("plot_template()",rcsid,__FILE__,__LINE__);
GR_report_error("Creating output file \"%s\"\n",filename);
GR_end_error();

fclose(fps);
return;
}


/* 
   This function makes nice plots showing the region of allowed masses
   in the tau0,tau1 plane and the ellipses representing filters on that
   plane.  It plots a single page, and is called by a driving routine.
*/
void plot_template_postscript(FILE *fps,struct Scope Grid,int clip,
            float mag,int page,int number,char *filename,char *date) {
int nfilt,i;
float mmn,mmx,mass;
float theta,dx0,dx1,pf,scalex,scaley,angle;
double tau0,tau1;
float maxrad,minrad,tau0min,tau0max,tau1min,tau1max;

/* find the number of templates */
nfilt=Grid.n_tmplt;

/* minimum and maximum mass of component stars */
mmx = Grid.m_mx ;
mmn = Grid.m_mn ;

/* angle of template ellipses relative to tau0,tau1 coordinates */
theta = Grid.theta;

/* minor/major diameter of template ellipses */
dx0 = Grid.dp ;
dx1 = Grid.dq ;

/* pi times the frequency */
pf=M_PI*(Grid.f_start) ;

/* find the maximum radius of the filters */
maxrad=minrad=Grid.dp;
if (Grid.dq>maxrad) maxrad=Grid.dq;
if (Grid.dq<minrad) minrad=Grid.dq;

/* step through templates to determine minimum/maximum tau values */
tau0min=tau1min=1.e30;
tau0max=tau1max=-1.e30;

for (i=0;i<nfilt;i++) {
	tau0=Grid.templates[i].tau0;
	tau1=Grid.templates[i].tau1;

	if (tau0<tau0min) tau0min=tau0;
	if (tau1<tau1min) tau1min=tau1;
	if (tau0>tau0max) tau0max=tau0;
	if (tau1>tau1max) tau1max=tau1;
}

/* rotation angle which make box "horizontal" */
angle=atan2(tau1max-tau1min,tau0max-tau0min);

/* scaling factors */
scalex=(PS_PW)/(tau1max-tau1min);
scaley=(PS_PW)/(tau0max-tau0min);

/* rotate the ouptput by this angle so that it does not angle up
   off the page. */
angle=(180.0/M_PI)*atan((tau1max-tau1min)/(tau0max-tau0min));

/* output postscript header into file */
fprintf(fps,"%%%%Page: graph %d\n",page);
fprintf(fps,"%% This is the top of page %d\n",page);
fprintf(fps,"%% There are several numbers below that can be edited by you to magnify, shrink or\n");
fprintf(fps,"%% shift the region being plotted.\n");
if (clip) {
	fprintf(fps,"/cpath {newpath %f %f moveto\n",PS_LM,PS_BM);
	fprintf(fps,"%f %f lineto\n",PS_LM,PS_BM+PS_CH);
	fprintf(fps,"%f %f lineto\n",PS_LM+PS_PW,PS_BM+PS_CH);
	fprintf(fps,"%f %f lineto\n",PS_LM+PS_PW,PS_BM);
	fprintf(fps,"%f %f lineto\n",PS_LM,PS_BM);
	fprintf(fps,"closepath} def\n");
	fprintf(fps,"%% COMMENT OUT THE FOLLOWING LINE TO DISABLE CLIPPING ON THIS PAGE\n");
	fprintf(fps,"cpath clip\n");
}
fprintf(fps,"/Times-Roman findfont 12 scalefont setfont\n");
fprintf(fps,"36 700 moveto (Template plot (file %s, page %d)) show\n",filename,page);
fprintf(fps,"36 685 moveto (Templates number %d) show\n",nfilt);
fprintf(fps,"36 670 moveto (%s) show\n",date);
fprintf(fps,"%f %f translate\n",PS_OX-(page-1)*PS_PW,PS_OY);
fprintf(fps,"%f rotate\n",-90-angle);
fprintf(fps,"%% UNCOMMENT THE FOLLOWING LINE TO FILL PAGE WITH FIGURE, OR ALTERNATIVELY\n");
fprintf(fps,"%%%f %f scale\n",-scalex,scaley);

/* choose smaller scale to preserve aspect ratio and still fit on page */
if (scalex<scaley) scaley=scalex; else scalex=scaley;

/* continue postscript header */
fprintf(fps,"%%...UNCOMMENT THE FOLLOWING LINE FOR 1:1 ASPECT RATIO.\n");
fprintf(fps,"%f %f scale\n",-scalex,scaley);
fprintf(fps,"1 setlinewidth\n");
fprintf(fps,"/el {gsave translate %f rotate %f %f scale\n",-180*theta/M_PI,0.5*dx1,0.5*dx0);
fprintf(fps,"                      0 0 1 0 360 arc stroke grestore} def\n");
fprintf(fps,"/num {gsave translate %f rotate 0 0 moveto\n",-180*theta/M_PI);
fprintf(fps,"/Times-Roman findfont %f scalefont setfont -1 1 scale dup\n",minrad);
fprintf(fps,"stringwidth pop -0.5 mul %f rmoveto show stroke -1 1 scale\n",-0.3*minrad);
fprintf(fps,"grestore } def\n");
fprintf(fps,"/ell {2 copy el num} def\n");
fprintf(fps,"%% REPLACE THE FOLLOWING NUMBER TO RESCALE (>1 MAGNIFIES)\n");
fprintf(fps,"%f\ndup scale\n",mag);
fprintf(fps,"%% REPLACE THE FOLLOWING TWO NUMBERS TO SHIFT IMAGE BY (tau0,tau1)\n");
fprintf(fps,"0 0\n-1 mul exch -1 mul translate\n");
fprintf(fps,"0.000001 setlinewidth\n");

/* draw a small set of tau0,tau1 axes */
fprintf(fps,"%% Starting small set of tau0, tau1 axes.\n");
fprintf(fps,"%f %f moveto\n",0.5*maxrad,0.5*maxrad);
fprintf(fps,"%f %f lineto stroke\n",1.3*maxrad,0.5*maxrad);
fprintf(fps,"/Symbol findfont %f scalefont setfont \n",0.2*maxrad);
fprintf(fps,"%f %f moveto 1 -1 scale (t1) show 1 -1 scale\n",1.5*maxrad,0.5*maxrad);

fprintf(fps,"%f %f moveto\n",0.5*maxrad,0.5*maxrad);
fprintf(fps,"%f %f lineto stroke\n",0.5*maxrad,1.3*maxrad);
fprintf(fps,"/Symbol findfont %f scalefont setfont \n",0.2*maxrad);
fprintf(fps,"%f %f moveto 1 -1 scale -90 rotate (t0) show 90 rotate 1 -1 scale\n",0.5*maxrad,1.5*maxrad);

/* draw graph outlining region satisfying mass constraints */
/* start point */
tau_of_mass(mmn,mmn,pf,&tau0,&tau1);

fprintf(fps,"%% Starting outline of mass-constraint region\n");
fprintf(fps,"newpath %f %f moveto\n",tau1-tau1min,tau0-tau0min);


/* equal mass curve */
fprintf(fps,"%% Equal mass curve (m1=m2=min to m1=m2=max)\n");
for (mass=mmn;mass<=mmx;mass+=0.001*(mmx-mmn)) {
	tau_of_mass(mass,mass,pf,&tau0,&tau1);
	fprintf(fps,"%f %f lineto\n",tau1-tau1min,tau0-tau0min);
}

/* maximum total mass curve */
fprintf(fps,"%% Maximum mass curve (m1=m2=max to m1=min,m2=max)\n");
for (mass=mmx;mass>=mmn;mass-=0.001*(mmx-mmn)) {
	tau_of_mass(mass,mmx,pf,&tau0,&tau1);
	fprintf(fps,"%f %f lineto\n",tau1-tau1min,tau0-tau0min);
}

/* minimum mass curve */
fprintf(fps,"%% Minimum mass curve (m1=min,m2=max to m1=min,m2=min)\n");
for (mass=mmx;mass>=mmn;mass-=0.001*(mmx-mmn)) {
	tau_of_mass(mass,mmn,pf,&tau0,&tau1);
	fprintf(fps,"%f %f lineto\n",tau1-tau1min,tau0-tau0min);
}

/* and draw the boundary of the mass constraint region */
fprintf(fps,"closepath stroke\n");

/* now step through the filters */
for (i=0;i<nfilt;i++) {
	tau0=Grid.templates[i].tau0;
	tau1=Grid.templates[i].tau1;
	if (number)
		fprintf(fps,"(%d) %f %f ell\n",i,tau1-tau1min,tau0-tau0min);
	else
		fprintf(fps,"%f %f el\n",tau1-tau1min,tau0-tau0min);

}

/* finish the postscript output */
fprintf(fps,"showpage\n");

return;	
}
