/* GRASP: Copyright 1997,1998  Bruce Allen */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "FrameL.h"
#include "grasp.h"
#define OLDNAMES 0            /* set to zero to use new channel names, 1 for old names */
#define LOCKLO 1
#define LOCKHI 10
#define CORRECTTIMESTAMPS 1   /* set to 1 to correct loss of timestamp resolution */
/* The compression method: None = 0; GZIP = 1; Diff = 2; Diff+GZIP = 3 */
#define COMPRESSION 3
/* The Level of GZIP compression used; Values between 1 and 9 allowed */
#define GZIP_LEVEL 1
/* the latest version of the frame lib that this code has been tested with */
#define FRAMELIB_TESTED 3.85

/* Each block of old-format data contains 5.07 secs of data.  This
parameter determines how many of these old-format blocks (now a frame)
end up in each FRAME file. */
#define FRAMES_PER_FILE 32
/* earth's equatorial radius, in meters */
#define EQUATORIAL (6.37814e+06)
/* earth's ellipticity or flattening due to rotation */
#define FLAT (3.35281e-3) 

/* the conversion from ADC counts to volts is: */
static char units[]="Units are 10 volts per 2048 counts.  Range -2048 to +2047";

#if (OLDNAMES)
/* channel assignments before Nov 15, 1994 */
static char *prenov15[]={
   "IFO output", "", "", "microphone", "dc strain", "mode cleaner pzt",
   "seismometer", "", "", "", "TTL locked", "arm 1 visibility", "arm 2 visibility",
   "mode cleaner visibility", "slow pzt", "arm 1 coil driver"};

/* channel assignments after Nov 15, 1994 */
static char *postnov15[]={
   "IFO output", "magnetometer", "microphone", "", "dc strain", "mode cleaner pzt",
   "seismometer", "slow pzt", "power stabilizer", "",
   "TTL locked", "arm 1 visibility", "arm 2 visibility", "mode cleaner visibility",
   "", "arm 1 coil driver"};
#else
/* channel assignments before Nov 15, 1994 */
static char *prenov15[]={
   "IFO_DMRO", "", "", "IFO_Mike", "IFO_DCDM", "PSL_MC_V",
   "IFO_Seis_1", "", "", "", "IFO_Lock", "IFO_EAT", "IFO_SAT",
   "IFO_MCR", "IFO_SPZT", "SUS_EE_Coil_V"};

/* channel assignments after Nov 15, 1994 */
static char *postnov15[]={
   "IFO_DMRO", "IFO_Mag_x", "IFO_Mike", "", "IFO_DCDM", "PSL_MC_V",
   "IFO_Seis_1", "IFO_SPZT", "PSL_PSS", "",
   "IFO_Lock", "IFO_EAT", "IFO_SAT", "IFO_MCR",
   "", "SUS_EE_Coil_V"};
#endif

/* Program's only argument is the name of the directory containing old-format data */
int main(int argc,char* argv[]) {
   char filename[256],name[256],hist[1024],*histnew,*buff,**chan_name;
   static char machinename[256]="";
   int i,code=1,num,large=50000,small=5000,n,first=1,firsttime=0,nlines;
   long buffSize;
   float fastrate=9868.4208984375,tblock,slowrate=986.84208984375,*real,*imag,*freq;
   double firstmsec=0.0,first_estimate,second_estimate,diff,dt,dtslow;
   float starttime=-100.0,guesstime;
   double currenttime=-200;
   int blockcount=0,channelsopen=0,expected;
   struct FrFile *outputfile=NULL;
   struct FrameH *frame;
   struct FrAdcData *adc[16];
   struct FrDetector *frdetect;
   struct FrVect *framevec;
   struct FrVect *framevecS;
   struct FrStatData *staticdata;
   struct FrStatData *staticdataS;
   struct ld_binheader bin_header;
   struct ld_mainheader main_header;
   struct tm timetm,*gtime,gts;
   time_t translate_time,calendartime;
   FILE *fp[16],*fpsweptsine,*pipe;
   void unhappyexit(int i);
   int get_run_number(int firsttime);

   /* print out some information about the library being used */
   fprintf(stderr,"translate compiled with Frame header file version: FRAMELIB_VERSION=%.2f\n",
	FRAMELIB_VERSION);
#if (FRAMELIB_VERSION_INT>=370)
   fprintf(stderr,"translate linked with Frame library archive version: FrLibVersion=%.2f\n",
	FrLibVersion(NULL));
   if ((int)(1000*(FRAMELIB_VERSION-FrLibVersion(NULL))))
      fprintf(stderr,
	"WARNING: translate code linked to different run-time library than header file version!\n");
#endif
   if (FRAMELIB_VERSION_INT!=100*FRAMELIB_VERSION)
      fprintf(stderr,
	"WARNING: in building this code FRAMELIB_VERSION_INT=%d != 100 x (FRAMELIB_VERSION=%.2f)\n",
         FRAMELIB_VERSION_INT,FRAMELIB_VERSION);
   if (FRAMELIB_VERSION>FRAMELIB_TESTED)
         fprintf(stderr,"Warning: translate has only been tested with FRAMELIB_VERSION <= %.2f\n",
	FRAMELIB_TESTED);

   /* initialize the frame system */
   FrLibIni(NULL,NULL,2);  
   buffSize=1000000;
   buff=malloc(buffSize);

   /* create a frame */
   frame=FrameHNew("C1");

   /* assign detector structure: site location and orientation information */
#if (FRAMELIB_VERSION_INT<=237)
   frame->detectRec=FrDetectorNew("real");
   frdetect=frame->detectRec;
   frdetect->latitude=34.1667;
   frdetect->longitude=118.133;
   frdetect->arm1Angle=180.0;
   frdetect->arm2Angle=270.0;
   frdetect->arm1Length=38.5;
   frdetect->arm2Length=38.1;
#else
   frame->detectProc=FrDetectorNew("real");
   frdetect=frame->detectProc;
   frdetect->latitudeD=34;
   frdetect->latitudeM=10;
   frdetect->latitudeS=0;
   frdetect->longitudeD=118;
   frdetect->longitudeM=8;
   frdetect->longitudeS=0;
   frdetect->armXazimuth=180.0;
   frdetect->armYazimuth=270.0;
   frdetect->armLength=38.3;
#endif

   /* Correct for oblateness of earth, use reference spheriod with
   flattening FLAT; EQUATORIAL is earth equatorial radius in meters.
   Reference: eqns (4.13-14) in "Spacecraft attitude determination and
   control", Ed. James R. Wortz, D. Reidel Publishing Co., Boston, 1985.
   Note: this SHOULD be corrected to add in the height of Caltech above
   sea level. */
   /* angle measured down from the North pole */

#if (FRAMELIB_VERSION_INT<=237)
   {
       float theta;
       theta=(M_PI/180.0)*(90.0-frdetect->latitude);
       frdetect->altitude=EQUATORIAL*(1.0-FLAT*cos(theta)*cos(theta));
   }
#else
       frdetect->elevation=0.0 /*FILL IN THE CORRECT VALUE */;
#endif

   /* now open files containing 40 meter data */
   if (!argv[1] || argc!=2) unhappyexit(1);

   /* step through all possible channels, seeing which channels have data */
   for (i=0;i<16;i++) {
      sprintf(name,"%s/channel.%d",argv[1],i);
      fp[i]=fopen(name,"r");
      if (fp[i]==NULL)
         fprintf(stderr,"File %s unavailable.  Skipping it...\n",name);
      else
         channelsopen++;
   }

   /* if there are no open files, then please exit with a warning message */
   if (channelsopen==0) unhappyexit(1);

   /* the sample times for the fast/slow channels */
   dt=1.0/fastrate;
   dtslow=1.0/slowrate;

   /* Define 4 fast, 12 slow ADC channels (long strings of blanks needed - see below) */
   for (i=0;i<16;i++)
      if (fp[i]!=NULL)
         if (i<4)
            /* sample rates differ from fastrate, slowrate -- see GRASP manual for details */
            adc[i]=FrAdcDataNew(frame,"                            ",50000.0*15.0/76.0,large,16);
          else
            adc[i]=FrAdcDataNew(frame,"                            ",5000.0*15.0/76.0,small,16);

   /* now loop over the input data, creating blocks of output data */
   while (code>0) {
      /* read a block of data */
      for (i=0;i<16;i++) {
         /* set size of data block */
         n=(i<4)?large:small;
         /* read data into frame short array */
         if (i<4 && fp[i]!=NULL)
            code=read_block(fp[i],&(adc[i]->data->dataS),&num,&tblock,&fastrate,0,&n,0,
                            &bin_header,&main_header);
         else if (fp[i]!=NULL)
            code=read_block(fp[i],&(adc[i]->data->dataS),&num,&tblock,&slowrate,0,&n,0,
                            &bin_header,&main_header);
      }

      /* if no data remains, we have found an error */
      if (code==0) {
         fprintf(stderr,"Error in translation: unexpected end of data!\n");
         abort();
      }

      /* check the various sample times */
      if (dt!=1.0/fastrate) fprintf(stderr,"Fast sample rates don't match!\n");
      if (dtslow!=1.0/slowrate) fprintf(stderr,"Slow sample rates don't match!\n");

      /* set time stamps for this block of data */
      /* Put the local California time-of-day into a structure for later use */
      timetm.tm_sec=main_header.tod_second;
      timetm.tm_min=main_header.tod_minute;
      timetm.tm_hour=main_header.tod_hour;
      timetm.tm_mday=main_header.date_day;
      timetm.tm_mon=main_header.date_month;
      timetm.tm_year=main_header.date_year;
      timetm.tm_wday=main_header.date_dow;
      timetm.tm_yday=-1;  /* info not available, but filled in by mktime */
      timetm.tm_isdst=-1; /* info not available, but filled in by mktime */

      /* Now convert the stored Calendar time into the right data type */
      calendartime=main_header.epoch_time_sec;

      /* Put the UTC time-of-day into a structure for later use */ 
      gtime=gmtime(&calendartime);
      gts=*gtime;

      /* set the time stamp for the first data sample (more precise than header time) */
      if (first) {
         firsttime=main_header.epoch_time_sec;
         firstmsec=0.001*main_header.epoch_time_msec;
         printf("UTC (gmtime) start time: %s",asctime(&gts));
         printf("          CA start time: %s\n",asctime(&timetm));

         /* assign the run number from 1,..,11 to the frame. */
         frame->run=get_run_number(firsttime);
         if (frame->run<1 || frame->run>11) unhappyexit(2);

         /* assign proper name to adc channel (to overwrite long blank space above) */
         if (frame->run<=2) {
            chan_name=prenov15;
            expected=11;
         }
         else {
            chan_name=postnov15;
            expected=13;
         }

         if (channelsopen!=expected) {
            fprintf(stderr,"Only found %d channels.  Expected %d\n",channelsopen,expected);
            exit(1);
         }

         for (i=0;i<16;i++)
            if (fp[i]!=NULL) {
               /* verify that name is correct */
               if (strcmp(chan_name[i],"")==0) {
                  fprintf(stderr,"Channel %d is not recognized and has no name!\n",i);
                  exit(1);
               }

               /* point to the correct channel name for this particular date, channel */
               strcpy(adc[i]->name,chan_name[i]);

               /* put in the physical volts/counts conversion */
#if (FRAMELIB_VERSION_INT<=237)
               adc[i]->data->unit[0]=(char *)malloc((strlen(units)+1)*sizeof(char));
               strcpy(adc[i]->data->unit[0],units);
#else
               adc[i]->data->unitY=(char *)malloc((strlen(units)+1)*sizeof(char));
               strcpy(adc[i]->data->unitY,units);
#endif

               /* which ADC "crate" was this */
               adc[i]->crate=i;
            }
      }

      if (CORRECTTIMESTAMPS) {
         guesstime=currenttime+76.0/15.0;
         if (fabs(guesstime-tblock)>1.0) {
            starttime=tblock;
            blockcount=0;
         }
         currenttime=(blockcount++)*((double)76.0/15.0)+starttime;

         /* put the time stamp into the frame structure */
         currenttime+=firstmsec;
#if (FRAMELIB_VERSION_INT<330)
         frame->UTimeS=firsttime+(int)currenttime;
         frame->UTimeN=(int)(1.e9*(currenttime-(int)currenttime));
         frame->dt=76.0/15.0;
#else
         frame->GTimeS=firsttime+(int)currenttime-UTCTOGPS;
         frame->GTimeN=(int)(1.e9*(currenttime-(int)currenttime));
         		                                /*JKB: should be INT(TAI-UTC) */
         frame->ULeapS=29;                              /* BA -- for Nov 1994 -- see GRASP manual on time defs */
         frame->dt=(76.0/15.0);
#endif

      }
      else {
         /* put the time stamp into the frame structure */
         tblock+=firstmsec;
#if (FRAMELIB_VERSION_INT<330)
         frame->UTimeS=firsttime+(int)tblock;
         frame->UTimeN=(int)(1.e9*(tblock-(int)tblock));
         frame->dt=num/slowrate;
#else
         frame->GTimeS=firsttime+(int)tblock-UTCTOGPS;
         frame->GTimeN=(int)(1.e9*(tblock-(int)tblock));
         		                                /*JKB: should be INT(TAI-UTC) */
         frame->ULeapS=29;                              /* BA -- for Nov 1995 -- see GRASP manual on time defs */
         frame->dt=(num/slowrate);
#endif

      }

      /* Localtime - UTC time in seconds */
      frame->localTime=-8*3600;

      /* frame->type[0]=0; */

      /* put in the history information (only once per translation) */
      if (first) {
         first=0;
         histnew=hist;
         time(&translate_time);

         /* get the name of the local machine */
         pipe=popen("hostname","r");
	 if (pipe==NULL) {
	    /* if we can't open the pipe, then list machine name as unknown */
           strcpy(machinename,"hostname undetermined");
         }
         else
           fscanf(pipe,"%s",machinename);

         histnew+=sprintf(histnew,"\nTranslation carried out by:\n");
         histnew+=sprintf(histnew,"     login: %s\n",getenv("LOGNAME"));
         histnew+=sprintf(histnew,"     user: %s\n",getenv("USER"));
         histnew+=sprintf(histnew,"     translation machine name: %s\n",machinename);
         histnew+=sprintf(histnew,"     directory: %s\n",getenv("PWD"));
         histnew+=sprintf(histnew,"     datapath: %s\n",argv[1]);
         histnew+=sprintf(histnew,"     translation program name: %s\n",argv[0]);
         histnew+=sprintf(histnew,"     source code name: %s\n","translate.c");
         histnew+=sprintf(histnew,"     Frame library header file (compile) version: %.2f\n",FRAMELIB_VERSION);
#if (FRAMELIB_VERSION_INT>=370)
         histnew+=sprintf(histnew,"     Frame library archive (link) version: %.2f\n",FrLibVersion(NULL));
#endif
         histnew+=sprintf(histnew,"     translation date: %s\n",ctime(&translate_time));
         FrHistoryAdd(frame,hist);

         /* read the swept sine calibration files (only once per run) */
         sprintf(name,"%s/swept-sine.ascii",argv[1]);
         fpsweptsine=fopen(name,"r");
         read_sweptsine(fpsweptsine,&nlines,&freq,&real,&imag);

         /* copy swept sine calibration data into vector; see below for packing style */
#if (FRAMELIB_VERSION_INT<=237)
         framevec=FrVectNew(FR_VECT_F,1,3*nlines,1.0,"Vifo/Vcoil");
#else
         framevec=FrVectNew(FR_VECT_4R,1,3*nlines,1.0,"Vifo/Vcoil");
#endif
         for (i=0;i<nlines;i++) {
            framevec->dataF[i]=freq[i];
            framevec->dataF[i+nlines]=real[i];
            framevec->dataF[i+2*nlines]=imag[i];
         }

         /* then link the calibration data into the history structure */
#if (FRAMELIB_VERSION_INT<330)
         staticdata=FrStatDataNew("sweptsine",
             "swept sine calibration:\npacking: freq[i], real[i], imaginary[i]",
             frame->UTimeS,INT_MAX,1,framevec);
#else
         staticdata=FrStatDataNew("sweptsine",
             "swept sine calibration:\npacking: freq[i], real[i], imaginary[i]",
             frame->GTimeS,INT_MAX,1,framevec);
#endif

#if (FRAMELIB_VERSION_INT<=237)
         FrStatDataAdd(&frame->detectRec->sData,staticdata);
#elif (FRAMELIB_VERSION_INT<=330)
         FrStatDataAdd(&frame->detectProc->sData,staticdata);
#else
         FrStatDataAdd(frame->detectProc,staticdata);
#endif

        /* put in lock range (INCLUSIVE low->high) Rolf: if 0=unlock and 1=lock 
           then you need LOCKLO=LOCKHI=1
         */

#if (FRAMELIB_VERSION_INT<=237)
         framevecS=FrVectNew(FR_VECT_S,1,2,1.0,"adcCounts");
#else
         framevecS=FrVectNew(FR_VECT_2S,1,2,1.0,"adcCounts");
#endif

         framevecS->dataS[0]=LOCKLO;  /* smallest value at which we are still in lock */
         framevecS->dataS[1]=LOCKHI;  /* largest value at which we are still in lock */

         /* then link the lockrange data into the history structure */
#if (FRAMELIB_VERSION_INT<330)
         staticdataS=FrStatDataNew("locklo/lockhi",
             "lock range:\npacking: array[0]=locklo array[1]=lockhi",
             frame->UTimeS,INT_MAX,1,framevecS);
#else
         staticdataS=FrStatDataNew("locklo/lockhi",
             "lock range:\npacking: array[0]=locklo array[1]=lockhi",
             frame->GTimeS,INT_MAX,1,framevecS);
#endif

#if (FRAMELIB_VERSION_INT<=237)
         FrStatDataAdd(&frame->detectRec->sData,staticdataS);
#elif (FRAMELIB_VERSION_INT<=330)
         FrStatDataAdd(&frame->detectProc->sData,staticdata);
#else
         FrStatDataAdd(frame->detectProc,staticdataS);
#endif
      }

      /* is the time stamp for this data block consistent with start time+offset? */
#if (FRAMELIB_VERSION_INT<330)
      first_estimate=frame->UTimeS+1.e-9*frame->UTimeN;
#else
      first_estimate=frame->GTimeS+UTCTOGPS+1.e-9*frame->GTimeN;
#endif

      second_estimate=main_header.epoch_time_sec+1.e-3*main_header.epoch_time_msec;
      diff=first_estimate-second_estimate;
      if (fabs(diff)>0.002)
         fprintf(stderr,"Time stamps have drifted by %f msec!\n",diff);

      /* Increment frame counter (set to 1 for first frame of each run) */
      frame->frame++;

      /* Open Frame file (one file per FRAMES_PER_FILE frames) */
      if ((frame->frame%FRAMES_PER_FILE)==1) {
         /* set file name.  Note than month=1 to 12 not 0 to 11! */
/* Obsolete as of Aug 1998 -- new file name is GPS time */
#if (FRAMELIB_VERSION_INT<330)
	sprintf(filename,"C1-%02d_%02d_%02d_%02d_%02d_%02d",gts.tm_year,gts.tm_mon+1,
		gts.tm_mday,gts.tm_hour,gts.tm_min,gts.tm_sec);
#else
	sprintf(filename,"C1-%d.F",frame->GTimeS);
#endif

         printf("Filename: %s\n",filename);
#if (FRAMELIB_VERSION_INT<330)
         outputfile=FrFileONew(filename, NO, buff, buffSize);
#else
         outputfile=FrFileONew(filename, COMPRESSION, buff, buffSize);
         if (GZIP_LEVEL>0) 
             {
                 printf("Building frames with compression gzip level = %d\n",GZIP_LEVEL);
             }
         FrFileOSetGzipLevel(outputfile,GZIP_LEVEL);
#endif
      }

      /* un-comment to print a short snippet of each Frame onto the screen */
      /* FrameDump(frame, stdout, 2); */

      /* Write frame to file, */
      FrameWrite(frame, outputfile);

      /* Close file if finished with FRAMES_PER_FILE or no remaining data */
      if ((frame->frame%FRAMES_PER_FILE)==0 || code==-1)
         FrFileOEnd(outputfile);
   }

   /* Free frame memory and return */                              
   FrameFree(frame);
   return(0);
}

/* this routine is called if something is wrong */
void unhappyexit(int i) {
switch (i) {
   case 1:
      fprintf(stderr,
      "Syntax: \ntranslate directory\nwhere channel.* files may be found in directory\n");
      exit(1);
   case 2:
      fprintf(stderr,
      "The UTC does not appear to lie in the range of any data set!\n");
      exit(1);
   default:
      abort();
}
return;
}

/* number of secs after Jan 1 1970 UTC at which Nov 1994 runs began */
static int stimes[]={784880277,784894763,785217574,785233119,785250938,785271063,
                     785288073,785315747,785333880,785351969,785368428,785388248};

/* This routine looks at the epoch time (sec) and returns the run number (1-11) */
int get_run_number(int firsttime) {
   int i;

   for (i=0;i<12;i++) 
      if (firsttime<stimes[i]) break;

   return i;
}
