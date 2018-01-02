/* GRASP: Copyright 1997,1998  Bruce Allen */

#define _XOPEN_SOURCE
#include "grasp.h"
#include "FrameL.h"
#include <string.h>
#undef _XOPEN_SOURCE

/* the latest version of the frame lib that this code has been tested with */
#define FRAMELIB_TESTED 3.85

static char *rcsid="$Id: frameinterface.c,v 1.36 2000/01/22 03:17:24 ballen Exp $\
\n$Name: RELEASE_1_9_8 $";

/* code will complain if the time-stamps (in secs) between frames don't
   agree as well as this */
#define RESOLUTION 3.2e-3

/* if frames contain UTC rather than GPS time, work properly anyway.  Set to 0 if you want to
use frames from >= 2004 */
#define WORKONBADFRAMES 1
static int imaxarg1,imaxarg2;
static int iminarg1,iminarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1)>(imaxarg2)?(imaxarg1):(imaxarg2))
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1)<(iminarg2)?(iminarg1):(iminarg2))

/* Possible name types for frame files - add more if needed. */
#ifdef _WIN32
/* should really be passed on compile line with -D! */
#define FRAMELIB_VERSION_INT 385 
static char *headnames[]={"C1-*.F","H-*.F","H-*.T","L-*.F","L-*.T","C1-*"};
#else
static char *headnames[]={"C1-*.F","H-*.F","H-*.T","L-*.F","L-*.T","C1-*[0-9]"};
#endif

/* This routine calls the function char *files() each time is needs a file with
   frames in it.  It then returns the first npoint points of data from the ADC labelled
   chname.  This is either in-lock or not-in-lock, as the user wishes.
*/
static char *datacodes="CSDFILfdsuilc";
/* The coding of the data type in (*adc[]->data).type, eg 0 corresponds to char, 
1 corresponds to short, 2 to double etc These codings are defined in FrameL.h :

enum FrVectTypes {FR_VECT_C,    vector of char           
                  FR_VECT_2S,   vector of short          
                  FR_VECT_8R,   vector of double         
                  FR_VECT_4R,   vector of float          
                  FR_VECT_4S,   vector of int            
                  FR_VECT_8S,   vector of long           
                  FR_VECT_C8,   vector of complex float  
                  FR_VECT_C16,  vector of complex double 
                  FR_VECT_S,    vector of string         
                  FR_VECT_2U,   vector of unsigned short 
                  FR_VECT_4U,   vector of unsigned int   
                  FR_VECT_8U,   vector of unsigned long  
                  FR_VECT_1U,   vector of unsigned char  
                  FR_VECT_END}; vector of unsigned char  */

int fget_ch(struct fgetoutput *fgetoutput,struct fgetinput *fgetinput) {
	static short *dats,*datl;			/* pointers to start of data, lock in current frame (set once!)*/
	static int ns,nl;				/* number of data, lock in current frame (set once!) */
	static int cs,cl;				/* current positions in current frane */
	static int ratio;				/* ratio of data/lock number (set once!) */
	static char *fname;				/* name of the file to open to get frames */
	static char *buff;				/* buffer in which data is stored by Frame library */
	static struct FrFile *iFile;			/* currently open file containing frames */
	static struct FrameH *frame=NULL;		/* currently open frame */
	static struct FrAdcData *adc0,*adcl;		/* adc's for IFO output and lock */
	static struct FrVect *frvects,*frvectl;		/* vector structures for IFO output and lock */
	static struct FrVect *frvectlohi;		/* vector structures for IFO lock low/hi values */
	static struct FrVect *calib;			/* vector structures for calibration info */
	static struct FrStatData *staticdata;		/* pointer to swept sine structure */
	static struct FrStatData *staticdataS;		/* pointer to slow channel structure */
	static struct FrStatData *staticdataL;		/* pointer to lock low/high values */
	static int databreak=1;				/* is there a break in the data? */
	static int fileopen=0;				/* is a file currently open? */
	static double firsttimestamp,lasttimestamp;	/* time of beginning/end of frame */
	static int lastframe,lastrun;			/* frame numbers of last frame & run */
	static long buffSize,debuglevel=0;		/* buffer size, debug level */
	static int entryno=0;				/* number of times entered */
	int npoint;					/* number of points requested on entry */
	int numneeded;					/* number of points needed to fill buffer */
	int imax,ncopy,brokenlock=0;			/* utility variables */
	int setlock=1;					/* if set to zero, indicates that we lost (and regained) lock */
	int top,i,start,end,where;			/* where to end the search for lock */
	static short **dat;				/* pointer to different data channels */
	static int *num,*curr,*rats;			/* number of points in different data channels,current pointer, ratios */
	static struct FrAdcData **adc;			/* adcs for other channels */
	static struct FrVect **frvect;			/* vector structures for other channels */
	static int nchanm1;				/* number of channels minus one */
	static int havestamped=0;			/* have we put the first timestamp on the data yet? */
	static short locklow,lockhi;			/* smallest value/largest value for which we are in lock */
	int fdesc;					/* file descriptor for the next file that we need to open */
        static int utctogps=UTCTOGPS;			/* store the offset between UTC and GPS times */
        static int timecheck=1;				/* check frame time-stamps for UTC vs GPS first time through */
#if (FRAMELIB_VERSION_INT>=340)
	static struct FrIO frIO;                        /* Frame3.40 uses structure to hold file descriptor */
#endif
	time_t printint;				/* used in printing the time */
	struct FrVect *slow_ch;                         /* vector structures for slow channel name info */


	/* this is the number of points of output that we need to generate */
	numneeded=npoint=fgetinput->npoint;

	/* count number of entries into this routine */
	if (++entryno==1) {
		/* if first entry, initialize the frame library, buffers */
		GR_start_error("fget_ch()",rcsid,__FILE__,__LINE__);
		GR_report_error("GRASP library linked with Frame header file version: FRAMELIB_VERSION=%.2f\n",FRAMELIB_VERSION);
#if (FRAMELIB_VERSION_INT>=370)
		GR_report_error("Executable linked with Frame library archive version: FrLibVersion=%.2f\n",FrLibVersion(NULL));
#endif
		if (FRAMELIB_VERSION_INT!=100*FRAMELIB_VERSION) {
			GR_report_error("WARNING: in building this code FRAMELIB_VERSION_INT=%d != 100 x (FRAMELIB_VERSION=%.2f)\n",
				FRAMELIB_VERSION_INT,FRAMELIB_VERSION);
			GR_report_error("This problem will arise if there is an error in the Makefile when\n");
			GR_report_error("GRASP is built, or if the FRAMELIB_VERSION in the FRAME library\n");
			GR_report_error("is not of the form X.yy where X=integer y=digit.\n");
		}
		if (FRAMELIB_VERSION>FRAMELIB_TESTED)
			GR_report_error("Warning: GRASP has only been tested with FRAMELIB_VERSION <= %.2f\n",FRAMELIB_TESTED);
		GR_end_error();

		FrLibIni(NULL,stderr,debuglevel);
		buffSize=1000000;
		buff=(char *)malloc((unsigned int)buffSize);
		fileopen=0;
		ns=nl=cs=cl=0;
		brokenlock=1;
		fgetoutput->discarded=0;
		fgetoutput->lastlock=0.0;
		fgetoutput->lastlock_gps=0.0-UTCTOGPS;
		fgetoutput->lostlock=0.0;
		fgetoutput->lostlock_gps=0.0-UTCTOGPS;
		nchanm1=fgetinput->nchan-1;
		for (i=0;i<nchanm1;i++) {
			num=(int*)malloc(sizeof(int)*nchanm1);
			curr=(int*)malloc(sizeof(int)*nchanm1);
			rats=(int*)malloc(sizeof(int)*nchanm1);
			adc=(struct FrAdcData**)malloc(sizeof(struct FrAdcData*)*nchanm1);
			frvect=(struct FrVect**)malloc(sizeof(struct FrVect*)*nchanm1);
			num[i]=curr[i]=rats[i]=0;
			adc[i]=(struct FrAdcData*)NULL;
			frvect[i]=(struct FrVect*)NULL;
		}
		if (fgetinput->nchan>1) dat=(short**)malloc(sizeof(short*)*nchanm1);
			
	}

	/* zero the number of points returned */
	for (i=0;i<fgetinput->nchan;i++)
		fgetoutput->npoint[i]=0;

	/* main function loop asks: do we still need to return points? */
	while (numneeded>0) {

		/* do we have a open file? */
		if (!fileopen) {
			/* see if the file() function is defined */
			if (fgetinput->files!=NULL) {
				/* if so, call files() function to get the names */
				fname=(*fgetinput->files)();
				fgetoutput->filename=fname;
				if (fname==NULL) {
					/* since no file to open, return 0 */
					GR_start_error("fget_ch()",rcsid,__FILE__,__LINE__);
					GR_report_error("Could not open NULL file name.\n");
					GR_report_error("Had %d points; still need %d points to satisfy request for %d samples.\n",
						npoint-numneeded,numneeded,fgetinput->npoint);
					GR_report_error("Discarding %d points remaining in the previous frame(s).\n",npoint-numneeded);
					GR_end_error();
					fgetoutput->discarded+=(npoint-numneeded);
					fgetoutput->returnval=0;
					return 0;
				}
				/* now open the frame file */
				iFile=FrFileINew(fname,buff,buffSize);
			}
			else {
				/* if it is not defined */
				fdesc=(*fgetinput->filedes)();
				if (fdesc<0) {
					/* since no file to open, return 0 */
					GR_start_error("fget_ch()",rcsid,__FILE__,__LINE__);
					GR_report_error("could not open file descriptor < 0\n");
					GR_report_error("had %d points; still need %d points...\n",npoint-numneeded,numneeded);
					GR_report_error("Discarding %d points remaining in the previous frame(s).\n",npoint-numneeded);
					GR_end_error();
					fgetoutput->discarded+=(npoint-numneeded);
					fgetoutput->returnval=0;
					return 0;
				}
				/* now open the frame file */
#if (FRAMELIB_VERSION_INT<=330)
				iFile=FrFileINewFd(fdesc,buff,buffSize);
#else
                                frIO.fd = fdesc;
				iFile=FrFileINewFd(&frIO,buff,buffSize);
#endif

			}
			if (iFile==NULL) {
				GR_start_error("fget_ch()",rcsid,__FILE__,__LINE__);
				GR_report_error("could not open file (%s)\n",fname);
				GR_end_error();
				fgetoutput->returnval=0;
				return 0;
			}
			fileopen=1;
		}

		/*  If less than zero points in current frame, something is badly wrong! */
		if (ns-cs<0) {
			GR_start_error("fget_ch()",rcsid,__FILE__,__LINE__);
			GR_report_error("less than zero points (%d) in frame!\n",ns-cs);
			GR_end_error();
			abort();
		}
		/* if there are points in the current frame, then */
		else if (ns-cs>0) {

			/* if we are curious about locked sections, look to see how long we are in lock */
			if (fgetinput->inlock) {
				/* find next point in lock stream covering needed range  */
				top=(cs+numneeded)/ratio;
				if (top%ratio!=0)
					top++;
				imax=IMIN(top,nl);

				/* then scan along to see if we drop out of lock during that interval */
				while (cl<imax && datl[cl]>=locklow && datl[cl]<=lockhi)
					cl++;
			}

			/* if we are curious about locked sections, and we are not in lock, look for return to lock */
			if (fgetinput->inlock && cl<imax) {

				/* if this is the first time that we lost lock, save the time */
				if (setlock) {
					setlock=0;
					fgetoutput->lostlock=firsttimestamp+ratio*(cl-1)/fgetoutput->srate;
					fgetoutput->lostlock_gps=fgetoutput->lostlock-UTCTOGPS;
				}

				/* Out of lock.  Look for return to lock. */
				imax=nl;
				while (cl<imax && (datl[cl]<locklow || datl[cl]>lockhi))
					cl++;
				fgetoutput->discarded+=(npoint-numneeded);
				fgetoutput->discarded+=(ratio*cl-cs);
				numneeded=npoint;
				cs=ratio*cl;
				brokenlock=1;

				/* zero the number of points returned */
				for (i=0;i<fgetinput->nchan;i++)
					fgetoutput->npoint[i]=0;
			}
			
			/* either don't care if in lock, or we are in lock, so copy some data! */
			else {	
				/* time to copy some points */
				ncopy=(ns-cs<numneeded)?ns-cs:numneeded;
				ncopy=IMIN(ns-cs,numneeded);

				/* if just starting the buffer, mark the start time, time into run */
				if (numneeded==npoint) {
					fgetoutput->tstart=firsttimestamp+cs/fgetoutput->srate;
					fgetoutput->tstart_gps=fgetoutput->tstart-UTCTOGPS;
					fgetoutput->dt=fgetoutput->tstart-fgetoutput->tfirst;

					/* if lost & regained lock, mark the time */
					if (!setlock && fgetinput->inlock) {
						fgetoutput->lastlock=fgetoutput->tstart;
						fgetoutput->lastlock_gps=fgetoutput->lastlock-UTCTOGPS;
					}
				}

				/* copy the points into the output buffer */
				if (!fgetinput->seek)
				  memcpy((void *)( fgetinput->locations[0]+(npoint-numneeded)*
						   ( ( ((*adc0->data).wSize) <= sizeof(short) ) ? 1 : ((*adc0->data).wSize)/sizeof(short) ) ), 
					 (const void *)(dats+cs),(size_t)(((*adc0->data).wSize)*ncopy));
				fgetoutput->npoint[0]+=ncopy;
					


				for (i=0;i<nchanm1;i++) {
					/* pointers to first, last locations.  note: INTEGER division! */
					start=cs/rats[i];
					if (cs%rats[i]!=0) start++;
					end=(cs+ncopy-1)/rats[i];
					if (start-end>1) {
						GR_start_error("fget_ch()",rcsid,__FILE__,__LINE__);
						GR_report_error("confused start/end pointers %d and %d for channel %s\n",start,end,
							fgetinput->chnames[i+1]);
						GR_end_error();
						abort();
					}
					where=fgetoutput->npoint[i+1];
					if (!fgetinput->seek)
						memcpy((void *)(fgetinput->locations[i+1]+where*
					( ( ((*adc[i]->data).wSize) <= sizeof(short) ) ? 1 : ((*adc[i]->data).wSize)/sizeof(short) ) ), 
						       (const void *)(dat[i]+start),(size_t)(((*(adc[i])->data).wSize)*(end-start+1)));
					fgetoutput->npoint[i+1]+=end-start+1;
				}

				/* adjust pointers to the data, and to the lock signal if needed */
				cs+=ncopy;
				numneeded-=ncopy;
				if (fgetinput->inlock) cl=cs/ratio;
			}
		}
		/* There are no points remaining in this frame... */
		else {
			/* try to read a new frame... */
			if (frame!=NULL)
				FrameFree(frame);
			frame=FrameRead(iFile);
			if (frame==NULL) {
				FrFileOEnd(iFile);
				fileopen=0;
				fname=NULL;
			}
			else {
				/* get the calibration information if it is needed */
				if (fgetinput->calibrate) {
#if (FRAMELIB_VERSION_INT<=237)
					staticdata=(struct FrStatData*)FrStatDataFind(frame->detectRec,"sweptsine");
#elif (FRAMELIB_VERSION_INT<=323)
					staticdata=(struct FrStatData*)FrStatDataFind(frame->detectProc,"sweptsine",frame->UTimeS);
#else
					staticdata=(struct FrStatData*)FrStatDataFind(frame->detectProc,"sweptsine",frame->GTimeS);
#endif
					if (staticdata==NULL) {
						GR_start_error("fget_ch()",rcsid,__FILE__,__LINE__);
						GR_report_error("Unable to locate \"sweptsine\" history structure in FRAME\n");
						GR_report_error("It appears that there is no calibration information in these frames!\n");
						GR_report_error("Sorry -- aborting...\n");
						GR_end_error();
						abort();
					}
					calib = (struct FrVect *) staticdata->data;
					fgetoutput->frinum=calib->nData;

/* Depending upon frame version, time stored in UTC or GPS time */
#if (FRAMELIB_VERSION_INT<=323)
					fgetoutput->tcalibrate=staticdata->timeStart;
					fgetoutput->tcalibrate_gps=fgetoutput->tcalibrate-UTCTOGPS;
#else
					fgetoutput->tcalibrate=staticdata->timeStart+utctogps;
					fgetoutput->tcalibrate_gps=fgetoutput->tcalibrate-UTCTOGPS;
#endif

					/* point to the calibration data */
					fgetoutput->fri=calib->dataF;
				}

			

				/* find a pointer to the first channel */
				adc0=FrAdcDataFind(frame,fgetinput->chnames[0]);
				/* Code added by Adrian Ottewill May 1999 to access Hanford Slow Channel data */
				for (i=0;i<fgetinput->nchan;i++) {
				  if (strcmp(fgetinput->chnames[i],"Slow")==0) {   
#if (FRAMELIB_VERSION_INT<=237) 
				    staticdataS=(struct FrStatData*)FrStatDataFind(frame->detectRec,"Slow_Ch_Name"); 
#elif (FRAMELIB_VERSION_INT<=323) 
				    staticdataS=(struct FrStatData*)FrStatDataFind(frame->detectProc,"Slow_Ch_Name",frame->UTimeS); 
#else 
				    staticdataS=(struct FrStatData*)FrStatDataFind(frame->detectProc,"Slow_Ch_Name",frame->GTimeS); 
#endif 			
				    slow_ch=(struct FrVect *) staticdataS->data; 
				    fgetoutput->slow_names=slow_ch->data; 
				  }  
				}

				if (adc0==NULL) {
					GR_start_error("fget_ch()",rcsid,__FILE__,__LINE__);
					GR_report_error("Channel \"%s\" does not exist in FRAME\n",fgetinput->chnames[0]);
					GR_report_error("Sorry -- aborting...\n");
					GR_end_error();
					abort();
				}

				/* get a pointer to the vector structure */
				frvects=adc0->data;

				/* find out how many points are in this frame */
				ns=frvects->nData;

				
				/* return the type of data of the principal channel if not short (the default) */
				if (frvects->type != 1) 
				  fgetinput->datatype[0]=datacodes[frvects->type];

				/* set the ratios of channel rates for the principal channel */
				fgetoutput->ratios[0]=1;

				for (i=0;i<nchanm1;i++) {
					/* find a pointer to the other channels */
					adc[i]=FrAdcDataFind(frame,fgetinput->chnames[i+1]);
					if (adc[i]==NULL) {
						GR_start_error("fget_ch()",rcsid,__FILE__,__LINE__);
						GR_report_error("Channel \"%s\" does not exist in FRAME\n",fgetinput->chnames[i+1]);
						GR_report_error("Sorry -- aborting...\n");
						GR_end_error();
						abort();
					}

					/* get a pointer to the vector structure */
					frvect[i]=adc[i]->data;

					/* find out how many points are in this frame, get pointers */
					num[i]=frvect[i]->nData;
					dat[i]=frvect[i]->dataS;

					/* set the ratios of channel rates for the non-principal channels */
					fgetoutput->ratios[i+1]=rats[i]=ns/num[i];

					
					/* find out the type of data in the non-principal channel if not short (the default) */
					if (frvect[i]->type != 1) 
					  fgetinput->datatype[i+1]=datacodes[frvect[i]->type];
				}

				if (fgetinput->inlock) {
					/* find the channel containing the lock signal */
					adcl=FrAdcDataFind(frame,"IFO_Lock");
					if (adcl==NULL) {
						GR_start_error("fget_ch()",rcsid,__FILE__,__LINE__);
						GR_report_error("Channel \"IFO_Lock\" does not exist in FRAME\n");
						GR_report_error("Sorry -- aborting...\n");
						GR_end_error();
						abort();
					}
					/* the lock data */
					frvectl=adcl->data;

					/* number of points in the lock data stream */
					nl=frvectl->nData;

					/* the ratio of principal/lock sample rates */
					ratio=ns/nl;

					/* pointer to the lock data, and set the current lock counter to zero */
					datl=frvectl->dataS;
					cl=0;

					/* now get the high/low values at which we are still in lock */
#if (FRAMELIB_VERSION_INT<=237)
					staticdataL=(struct FrStatData*)FrStatDataFind(frame->detectRec,"locklo/lockhi");
#elif (FRAMELIB_VERSION_INT<=323)
					staticdataL=(struct FrStatData*)FrStatDataFind(frame->detectProc,"locklo/lockhi",frame->UTimeS);
#else
					staticdataL=(struct FrStatData*)FrStatDataFind(frame->detectProc,"locklo/lockhi",frame->GTimeS);
#endif
					if (staticdataL==NULL) {
						GR_start_error("fget_ch()",rcsid,__FILE__,__LINE__);
						GR_report_error("Unable to locate \"locklo/lockhi\" history structure in FRAME\n");
						GR_report_error("It appears that there is no lock range information in these frames!\n");
						GR_report_error("Sorry -- aborting...\n");
						GR_end_error();
						abort();
					}
					/* pointer to the array containing the low/high values */
					frvectlohi = (struct FrVect*) staticdataL->data;

					/* record the low/high values internally, and return them to the user also */
					fgetoutput->locklow=locklow=frvectlohi->dataS[0];
					fgetoutput->lockhi=lockhi=frvectlohi->dataS[1];
				}
				cs=0;

				/* sample rate (FIRST ONE NOT ACCURATE ENOUGH!) */
#if (FRAMELIB_VERSION_INT<=237)
				fgetoutput->srate=adc0->frequency;
#else
				fgetoutput->srate=adc0->sampleRate;
#endif

				/* and get a pointer to the shorts */
				dats=frvects->dataS;

#if (FRAMELIB_VERSION_INT<=237)
				firsttimestamp=(double)(frame->UTimeS)+1.e-9*(double)(frame->UTimeN);
#elif (FRAMELIB_VERSION_INT<=323)
				firsttimestamp=(double)(frame->UTimeS)+1.e-9*(double)(frame->UTimeN);
#else
                                /* Depending upon frame version, time stored in UTC or GPS time */
				firsttimestamp=(double)(frame->GTimeS)+1.e-9*(double)(frame->GTimeN);

				/* check that GPS time is what is actually stored! */
				if (timecheck && firsttimestamp>7.539e8 && WORKONBADFRAMES) {
					GR_start_error("fget_ch()",rcsid,__FILE__,__LINE__);
					printint=(time_t)firsttimestamp+UTCTOGPS;
					GR_report_error("The time stamps from these frames show a year >= 2004, APPROXIMATELY UTC %s",
							asctime(gmtime(&printint)));
					GR_report_error("I am assuming that these frames contain 1994 data with Unix-C rather than the\n");
					GR_report_error("(correct!) GPS time stamps, and am correcting appropriately, so that\n");
					printint=(time_t)firsttimestamp;
					GR_report_error("the data dates from APPROXIMATELY UTC %s",asctime(gmtime(&printint)));
					printint=(time_t)firsttimestamp+UTCTOGPS;
					GR_report_error("If this message is erroneous, and you are REALLY processing data from %s",
							asctime(gmtime(&printint)));
					GR_report_error("please recompile this file with WORKONBADFRAMES defined to 0.\n");
					GR_report_error("Documentation on time-stamp problems in some distributed frames is at:\n");
					GR_report_error("http://www.lsc-group.phys.uwm.edu/~ballen/grasp-distribution/faqs.html\n");
					GR_end_error();
                                        utctogps=timecheck=0;
				}
				firsttimestamp+=utctogps;
#endif

				/* set initial time, offset time */
				if (havestamped==0) {
					fgetoutput->tfirst=firsttimestamp;
					fgetoutput->tfirst_gps=fgetoutput->tfirst-UTCTOGPS;
					havestamped=1;
				}

				/* on first entry into program, set the loss-of-lock time stamp properly */
				if (entryno==1 && fgetinput->inlock && setlock) {
					setlock=0;
					fgetoutput->lostlock=firsttimestamp;
					fgetoutput->lostlock_gps=fgetoutput->lostlock-UTCTOGPS;
				}

				if (!databreak && (fabs(firsttimestamp-lasttimestamp)>RESOLUTION)) {
					GR_start_error("fget_ch()",rcsid,__FILE__,__LINE__);
					GR_report_error("FRAMES NOT SEQUENTIAL\n");
					GR_report_error("run %3d frame %5d  ended  at Unix-C time:   %f sec\n",lastrun,lastframe,lasttimestamp);
					GR_report_error("                                GPS time:   %f sec\n",lasttimestamp-UTCTOGPS);
					GR_report_error("run %3d frame %5d started at Unix-C time: %f sec\n",(int)frame->run,(int)frame->frame,firsttimestamp);
					GR_report_error("                                  GPS time: %f sec\n",firsttimestamp-UTCTOGPS);
					GR_report_error("Time gap is %f sec\n",firsttimestamp-lasttimestamp);
					GR_report_error("Gap starts %f seconds into run; ends %f seconds into run.\n",
						lasttimestamp-fgetoutput->tfirst,firsttimestamp-fgetoutput->tfirst);
					GR_report_error("Discarding %d points remaining in the previous frame(s).\n",npoint-numneeded);
					GR_end_error();
					fgetoutput->discarded+=(npoint-numneeded);
					numneeded=npoint;
					/* zero the number of points returned */
    					for (i=0;i<fgetinput->nchan;i++)
						fgetoutput->npoint[i]=0;
					databreak=1;
					brokenlock=1;
					if (setlock && fgetinput->inlock) {
						setlock=0;
						fgetoutput->lostlock=lasttimestamp;
						fgetoutput->lostlock_gps=fgetoutput->lostlock-UTCTOGPS;
					}


				}
				else	
					databreak=0;
				lasttimestamp=firsttimestamp+frame->dt;

				lastframe=frame->frame;
				lastrun=frame->run;
			}
		}
	}
	/* since numneeded is zero we are done */
	if (brokenlock==1) {
		fgetoutput->returnval=1;
		return 1;
	}
	else {
		fgetoutput->returnval=2;
		return 2;
	}
}

/* 
   method=0 => rational function interpolation
   method=1 => polynomial interpolation
   order in either case is the order of interpolation
*/

void GRcalibrate(float *fri,int frinum,int num,float *complex,float srate,int method,int order) {
 	void hunt(float xx[], unsigned long n, float x, unsigned long *jlo);
 	void polint(float xa[], float ya[], int n, float x, float *y, float *dy);
	void ratint(float xa[], float ya[], int n, float x, float *y, float *dy);
	void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);
	void splint(float xa[], float ya[], float y2a[], int n, float x, float *y);

	float *ssf,*ssr,*ssi,freq,interp_real,interp_imag,*ssr2d=NULL,*ssi2d=NULL;
	float error_real,error_imag;
	int i,k,ssn;
	unsigned long jlo=1;

       /* read in the swept sine file */
       ssn=frinum/3;
       ssf=fri;
       ssr=ssf+ssn;
       ssi=ssr+ssn;

if (ssn<2) {
                GR_start_error("GRcalibrate()",rcsid,__FILE__,__LINE__);
                GR_report_error("Less than 2 lines (frequency values) of calibration data!\n");
                GR_report_error("Something is wrong with the data (in the frame?)... exiting.\n");
                GR_end_error();
		exit(1);
}

/* if using spline interpolation, please set up 2nd derivatives */
if (method==2) {
	ssr2d=(float*)malloc(ssn*sizeof(float));
	ssi2d=(float*)malloc(ssn*sizeof(float));
	spline(ssf-1,ssr-1,ssn,2.e30,2.e30,ssr2d-1);
	spline(ssf-1,ssi-1,ssn,2.e30,2.e30,ssi2d-1);
}

/* step through num/2+1 different (non-negative) frequencies */
for (i=0;i<=num/2;i++) {

	/* frequency at which we want response function */
	freq=(i*srate)/num;

	if (method==2) {
		/* use spline interpolation to estimate value */
		splint(ssf-1,ssr-1,ssr2d-1,ssn,freq,&interp_real);	
		splint(ssf-1,ssi-1,ssi2d-1,ssn,freq,&interp_imag);
	}	
	else {
		/* hunt for nearest entries in table (NR routine) -1 for zero offset */
		hunt(ssf-1,ssn,freq,&jlo);

		/* find where to interpolate in table, see NR in C section 3.4 */
		k=IMIN(IMAX(jlo-(order-1)/2,1),ssn+1-order);

		/* now get the interpolated value (-2 is zero offset, twice) */
		if (method==0) {
			/* rational function interpolation */
			ratint(ssf+k-2,ssr+k-2,order,freq,&interp_real,&error_real);
			ratint(ssf+k-2,ssi+k-2,order,freq,&interp_imag,&error_imag);
		}
		else if (method==1) {
			/* polynomial interpolation */
			polint(ssf+k-2,ssr+k-2,order,freq,&interp_real,&error_real);
			polint(ssf+k-2,ssi+k-2,order,freq,&interp_imag,&error_imag);
		}
		else {
			GR_start_error("calibrate()",rcsid,__FILE__,__LINE__);
			GR_report_error("unrecognized interpolation ");
			GR_report_error("method %d\n",method);
			GR_end_error();
			abort();
		}
	}

	/* store the interpolated values in arrays */
	complex[2*i]=interp_real;
	complex[2*i+1]=interp_imag;
}

/* if using spline interpolation, free the associated memory */
if (method==2) {
	free(ssr2d);
	free(ssi2d);
}

return;
}





void GRnormalize(	float *fri,	/* array to get frequency, real, imaginary parts.        */
			int frinum,	/* number of points in previous array (divisible by 3)   */
			int npoint,	/* number of points in time domain sample (power of 2)  */
			float srate,	/* sample rate in Hz of the gw signal stream            */
			float *response /* on return, contains complex response function        */
) {

	float freq,factor,f2,mod2,ssreal,ssimag;
	int j,ir,ii;
	
	/* get interpolated swept sine SS curve (store in response[] for now) */
	GRcalibrate(fri,frinum,npoint,response,srate,2,0);

	/* assemble all the factors -- square root bandwidth */
	factor=sqrt(9.35);

	/* Bob Spero's constant k */
	factor/=21399.0;

	/* volts IFO/ADC bit */
	factor*=(10.0/2048.0);

	/* from two time derivatives */
	factor/=(-4.0*M_PI*M_PI);
	
	/* set DC response to zero */
	response[0]=response[1]=0.0;

	/* set Nyquist response to zero */
	response[npoint]=response[npoint+1]=0.0;

	/* loop over all frequencies except DC (j=0) & Nyquist (j=npoint/2) */
	for (j=1;j<npoint/2;j++) {
		/* subscripts of real, imaginary parts */
		ir=2*j;
		ii=ir+1;

		/* frequency and frequency^2 */
		freq=srate*j/npoint;
		f2=freq*freq;

		/* real, imaginary parts of swept sine SS */
		ssreal=response[ir];
		ssimag=response[ii];

		/* squared modulus of swept sine SS */
		mod2=ssreal*ssreal+ssimag*ssimag;

		/* response function (one over complex conjugate of SS) */
		response[ir]=factor*ssreal/(f2*mod2);
		response[ii]=factor*ssimag/(f2*mod2);
	}
	return;
}

#define MAXFILES 16384
#define MAXLEN 128
#ifdef REALTIME
#include "cadef.h"
#include "ezca.h"
#endif

char *framefiles() {
	static int entry=0,nfiles=0,realtime;
	char longpath[256],*pathhead;
	const char *environment_variable="GRASP_FRAMEPATH",*environment_variable2="GRASP_REALTIME";
        static char filelist[MAXFILES][MAXLEN];
        int namelen,nametype;
	FILE *fp;
        /*length of array of name-expressions */
	const int numheads=sizeof(headnames)/sizeof(char*);
#ifdef REALTIME
	static char filename[256]="a";
	int value;
#endif
	if (entry++==0) {
		/* see if we want real time data */
		pathhead=getenv(environment_variable2);

		if (pathhead) {
			realtime=1;
			/* tell user that we will use real-time data */
			GR_start_error("framefiles()",rcsid,__FILE__,__LINE__);
			GR_report_error("environment variable %s is set.\n",environment_variable2);
			GR_report_error("getting real-time data... type\n");
			GR_report_error("unsetenv %s\nif this is not what you want!\n");
#ifndef REALTIME
			GR_report_error("Sorry: this code compiled without the real-time option REALTIME defined.\n");
			GR_report_error("it will require re-compiliation, after you modify SiteSpecific or Makefile.\n");
			GR_end_error();  /* fix made by Bruce Allen, Nov 10, 1997 */
			exit(1);
#endif
			GR_end_error();  /* fix made by Bruce Allen, Nov 10, 1997 */
#ifdef REALTIME
			ezcaSetMonitor("fb_file_string",ezcaString);
			/* GR_report_error("********************* Calling ezcaSetMonitor ********************\n"); */
#endif
		}
		else {
			realtime=0;

			/* get the environment variable giving the data path */
			pathhead=getenv(environment_variable);

			/* if this is null, then print out a warning to the user */
			if (!pathhead) {
				GR_start_error("framefiles()",rcsid,__FILE__,__LINE__);
				GR_report_error("the environment variable %s must be set.\n",environment_variable);
				GR_report_error("It should point to a directory containing frame files.\n");
				GR_report_error("It can be set with a command like:\nsetenv %s /this/is/the/path\n\n",environment_variable);
#ifdef REALTIME
				GR_report_error("As an alternative, this version of the code permits real-time analysis.\n");
				GR_report_error("to enable this option, set the environment variable %s with the command:\n",environment_variable2);
				GR_report_error("setenv %s\n\n",environment_variable2);
#endif
				GR_end_error();
				exit(1);
			}
			/* construct commands to list file names: loop until some files are located */
			nametype=-1;
			while (nfiles==0 && nametype<numheads-1) {
				nametype++;

				/* pointer to the file names */
#ifdef _WIN32
				sprintf(longpath,"dir %s\\%s /b",pathhead,headnames[nametype]);
				fp=_popen(longpath,"r");
#else
				sprintf(longpath,"ls %s/%s 2>/dev/null",pathhead,headnames[nametype]);
				fp=popen(longpath,"r");
#endif
				if (fp==NULL) {
					GR_start_error("framefiles()",rcsid,__FILE__,__LINE__);
					GR_report_error("unable to open a pipe: popen(%s,\"r\") failed!\n",longpath);
					GR_end_error();
					exit(1);
				}

				/* read file names */
				while (EOF!=fscanf(fp,"%s\n",filelist[nfiles])) {
					namelen=strlen(filelist[nfiles]);
					if (namelen>MAXLEN-1) {
						GR_start_error("framefiles()",rcsid,__FILE__,__LINE__);
						GR_report_error("file name:\n%s\nis too long (%d characters)! MAXLEN = %d\n",
									filelist[nfiles],namelen,MAXLEN);
						GR_report_error("please recompile framefiles() with MAXLEN >= %d\n",namelen+1);
						GR_end_error();
						exit(1);
					}
					nfiles++;
					if (nfiles>MAXFILES) {
						GR_start_error("framefiles()",rcsid,__FILE__,__LINE__);
						GR_report_error("more than MAXFILES=%d files %s found in directory %s.\n",
							MAXFILES,headnames[nametype],environment_variable);
						GR_report_error("you will need to recompile framefiles() with larger value of MAXFILES\n");
						GR_end_error();
						exit(1);
					}
				}
#ifdef _WIN32 
				_pclose(fp);
#else 
				pclose(fp);
#endif 
			}
		
			if (nfiles==0) {
				GR_start_error("framefiles()",rcsid,__FILE__,__LINE__);
				GR_report_error("No files of forms:\n");
				for (nametype=0;nametype<numheads;nametype++)
					GR_report_error("\t%s\n",headnames[nametype]);
				GR_report_error("found in directory %s specified by the environment variable %s.\n",pathhead,environment_variable);
				GR_report_error("This environment variable should point to a directory containing frame files,\n");
				GR_report_error("and can be set with a command like:\nsetenv %s /this/is/the/path\n",environment_variable);
				GR_end_error();
				exit(1);
			}

			GR_start_error("framefiles()",rcsid,__FILE__,__LINE__);
			GR_report_error("Environment variable %s set to: %s\n",environment_variable,pathhead);
			GR_report_error("using %d files of type %s from directory %s\n",nfiles,headnames[nametype],pathhead);
			GR_end_error();			
		}
	}

#ifdef REALTIME
	if (realtime) {
		while (1!=(value=ezcaNewMonitorValue("fb_file_string",ezcaString))) {
			/* sleep here! */
			/* GR_report_error(" you asked me for this: %d....\n",value); */
		};
		ezcaGet("fb_file_string",ezcaString,1,filename);
		/* GR_report_error("Call %d to framefiles got filename: %s (value=%d)\n",entry,filename,value); */
		return(filename);
	}
#endif
	if (entry>nfiles)
		return NULL;
	else
		return filelist[entry-1];
}

