#include "gl.h"
#include "glx.h"
#include "grasp.h"

Window make_rgb_db_window( Display *, unsigned int, unsigned int);

void plottf_bw(float **pic,int pdim)
{

    int i,j;
    float magn;
    static Display *dpy;
    static Window win;
    Window make_rgb_db_window( Display *,unsigned int, unsigned int);
    static int first=1;
    
    if(first){
        first=0;    
        dpy = XOpenDisplay(NULL);
        win = make_rgb_db_window( dpy, pdim, pdim);
        XMapWindow( dpy, win );
        glClearColor (0.0, 0.0, 0.0, 0.0);
        glClear(GL_COLOR_BUFFER_BIT);
        glMatrixMode (GL_PROJECTION);
        glLoadIdentity ();
        glOrtho(0.0, pdim*1.0, 0.0, pdim*1.0, 0.0, pdim*1.0);
    }
	glBegin(GL_POINTS);
    for(i=0;i<pdim;i++){
        for(j=0;j<pdim;j++){
            magn=*(pic[i]+j);
            glColor3f(magn,magn,magn);
            glVertex2i(i,j);
        }
    }
	glEnd();
    glXSwapBuffers( dpy, win );
    glFlush();
}





void plottf(float **pic,int pdim)
{

    int i,j;
    float magn;
    static Display *dpy;
    static Window win;
    Window make_rgb_db_window( Display *,unsigned int, unsigned int);
    float colr,colg,colb;
    static int first=1;
    
    if(first){
        first=0;    
        dpy = XOpenDisplay(NULL);
        win = make_rgb_db_window( dpy, pdim, pdim);
        XMapWindow( dpy, win );
        glClearColor (0.0, 0.0, 0.0, 0.0);
        glClear(GL_COLOR_BUFFER_BIT);
        glMatrixMode (GL_PROJECTION);
        glLoadIdentity ();
        glOrtho(0.0, pdim*1.0, 0.0, pdim*1.0, 0.0, pdim*1.0);
    }
	glBegin(GL_POINTS);
    for(i=0;i<pdim;i++){
        for(j=0;j<pdim;j++){
            magn=*(pic[i]+j);
            huetorgb(magn,&colr,&colg,&colb);
            glColor3f(colr,colg,colb);
            glVertex2i(i,j);
        }
    }
	glEnd();
    glXSwapBuffers( dpy, win );
    glFlush();
}


Window make_rgb_db_window( Display *dpy,
				  unsigned int width, unsigned int height )
{
   int attrib[] = { GLX_RGBA,
		    GLX_RED_SIZE, 1,
		    GLX_GREEN_SIZE, 1,
		    GLX_BLUE_SIZE, 1,
		    GLX_DOUBLEBUFFER,
		    None };
   int scrnum;
   XSetWindowAttributes attr;
   unsigned long mask;
   Window root;
   Window win;
   GLXContext ctx;
   XVisualInfo *visinfo;

   scrnum = DefaultScreen( dpy );
   root = RootWindow( dpy, scrnum );

   visinfo = glXChooseVisual( dpy, scrnum, attrib );
   if (!visinfo) {
      printf("Error: couldn't get an RGB, Double-buffered visual\n");
      exit(1);
   }

   /* window attributes */
   attr.background_pixel = 0;
   attr.border_pixel = 0;
   attr.colormap = XCreateColormap( dpy, root, visinfo->visual, AllocNone);
   attr.event_mask = StructureNotifyMask | ExposureMask;
   mask = CWBackPixel | CWBorderPixel | CWColormap | CWEventMask;

   win = XCreateWindow( dpy, root, 0, 0, width, height,
		        0, visinfo->depth, InputOutput,
		        visinfo->visual, mask, &attr );

   ctx = glXCreateContext( dpy, visinfo, NULL, True );

   glXMakeCurrent( dpy, win, ctx );

   return win;
}









