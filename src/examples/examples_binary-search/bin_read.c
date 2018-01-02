#include <stdio.h>
#include <string.h>
#include "grasp.h"


#define PBINS     20     /* used to convert from rsq stat to chisq stat */
#define RSQ     2.48     /* rsq threshold */

extern char *optarg;
extern int optind;

const char *rcsid="$Id";

struct MassPair {
  float one;
  float two;
};

struct MF_Header {
  char *text;
  int one;
  int num_templates;
  float f_low;
  float srate;
  struct MassPair *mass;
};

struct MF_Signal {
  float distance;
  float snr_max;
  float snr_max_pow_renorm;
  float snr_max_med_renorm;
  int offset;
  float phase;
  float stat;
  int num_crossing;
};

struct MF_Segment {
  double datastart;
  int is_gauss;
  struct MF_Signal *signal;
};

/* byteorder_same == 1 if this machine has same byte-ordering
   (big vs little endian) as the machine that took the data.
   byteorder_same == 0 if they have opposite byte-ordering.
*/
int byteorder_same=-1;


int main(int argc, char *argv[])
{
  FILE *fopen_segment(char *, int);
  int read_header(FILE *, struct MF_Header *);
  int read_segment(FILE *, struct MF_Segment *, int);
  int event(float, float, int, struct MF_Segment);
  extern int getopt (int argc, char *const *argv, const char *shortopts);
  void usage(char *program);

  double lasttime=-1,time=0;
  float snr_max=-1, rsq = RSQ;
  int i,j,seg=-1,total=0;
  int logfile=0;
  char *mfpath, logfname[256];
  FILE *fp,*fplog;
  struct MF_Header header={0};

  /* get path to read the data */
  mfpath = getenv("GRASP_MFPATH");
  if (mfpath==NULL) {
    GR_start_error("main",rcsid,__FILE__,__LINE__);
    GR_report_error("the environment variable GRASP_MFPATH must be set.\n");
    GR_end_error();
    return 1;
  }

  /* parse the optional arguments */
  while (1) {
    int c;
      
    /* call standard C library option parser */
    c = getopt(argc, argv, "hlr:");
    if (c == -1) break;
      
    switch (c) {

    case 'h': /* Print a simple usage message */
      usage(argv[0]);
      exit(0);

    case 'l': /* Follow a logfile */
      logfile = 1;
      sprintf(logfname,"%s/insert.log\0",mfpath);
      fplog = fopen(logfname,"r");
      if (fplog == NULL) {
	GR_start_error("main()",rcsid,__FILE__,__LINE__);
	GR_report_error("Unable to open injection logfile");
	GR_end_error();
	exit(1);
      }
      break;

    case 'r': /* Set rsq threshold */
      rsq = atof(optarg);
      break;

    default: /* If things fail for some reason */
      fprintf(stderr,"?? getopt returned character code 0%o ??", c);
    }
  }
  
  /* Deal with possible errors */
  if (optind < argc) {
    fprintf(stderr,"non-option ARGV-elements: ");
    while (optind < argc) fprintf(stderr,"%s ", argv[optind++]);
  }

  /* Read the header file */
  fp = fopen_segment(mfpath,-1);
  if (!read_header(fp,&header)) {
    GR_start_error("main()",rcsid,__FILE__,__LINE__);
    GR_report_error("Unable to read header file");
    GR_end_error();
    exit(1);
  }
  fclose(fp);


  while (1) {

    struct MF_Segment segment={0};
    struct MF_Signal *signal;

    /* get segment and time of next injection */
    if (logfile) {
      lasttime = time;
      if (EOF==fscanf(fplog,"%d %*d %lf %*f %*f %*f %*f\n",&seg,&time))
	break;
    } else {
      seg++;
    }

    if (seg % 100 == 0) fprintf(stderr," [%d]\r",seg);

    /* open the appropriate segment file */
    if (NULL==(fp=fopen_segment(mfpath,seg)))
      break;

    /* read data file */
    if (!read_segment(fp,&segment,header.num_templates)) {
      GR_start_error("main()",rcsid,__FILE__,__LINE__);
      GR_report_error("Unable to read segment %05d file",seg);
      GR_end_error();
      exit(1);
    }
    fclose(fp);

    /* if (!segment.is_gauss)
         continue; */

    if (time != lasttime) {
      ++total;
      if (snr_max >= 0) printf("%f\n",snr_max);
      snr_max = 0;
    }

    for (i = 0; i < header.num_templates; ++i)
      if (segment.signal[i].stat < rsq)
        if (segment.signal[i].snr_max > snr_max)
          snr_max = segment.signal[i].snr_max;

    signal = segment.signal;
    free(signal);
  }
  
  /* last one... */
  ++total;
  if (snr_max >= 0) printf("%f\n",snr_max);

  /* close logfile */
  if (logfile) fclose(fplog);

  fprintf(stderr,"total: %d\n",total);

  return 0;
}


/* Open file for the specified segment; open header if segment negative */
FILE *fopen_segment(char *path, int segment)
{
  static char fname[128];
  if (segment<0) sprintf(fname,"%s/signal.header",path);
  else sprintf(fname,"%s/signal.%05d",path,segment);
  return fopen(fname,"r");
}


/* This routine swaps byte order on length-n arrays of size-byte objects */
void swap(int size,int n,char *input) {
  char temp;
   int i,j,hi,lo,top,bot;

  /* for each element of the array */
  for (i=0;i<n;i++) {
    lo=i*size;
    hi=lo+size-1;
    /* copy the element in reverse order into storage */
    for (j=0;j<size/2;j++) {
      top=lo+j;
      bot=hi-j;
      temp=input[top];
      input[top]=input[bot];
      input[bot]=temp;
    }
  }
  return;
}


#define DEFAULT_HEADER_SIZE 4096
/* Read the header file --- see doumentation */
int read_header(FILE *fp, struct MF_Header *header)
{
  int bytes=0,headerlength,code;
  /* read the size of the header from the first characters in the file */
  code=fscanf(fp,"%d\n",&headerlength);
  /* if first characters don't contain the header length, use default size */
  if (code!=1) {
        headerlength=DEFAULT_HEADER_SIZE;
  }
  /* rewind file to the beginning */
  rewind(fp);
  header->text = (char *)realloc((void *)header->text,headerlength);
  if (header->text==NULL) {
    GR_start_error("read_header()",rcsid,__FILE__,__LINE__);
    GR_report_error("Unable to allocate %d bytes of memory\n",
                    headerlength);
    GR_end_error();
    return 0;
  }
  bytes += fread(header->text,1,headerlength,fp);
  bytes += fread(&header->one,1,4,fp);
  /* determine byte-ordering of machine */
  byteorder_same=(header->one==1)?1:0;
  /* check that we got it right */
  if (!byteorder_same) {
    swap(4,1,(char *)&header->one);
    if (header->one!=1) {
      GR_start_error("read_header()",rcsid,__FILE__,__LINE__);
      GR_report_error("Unable to determine byte-order of machine!\n");
      GR_report_error("It does not appear to be either"
                      "big- or little-endian!\n");
      GR_end_error();
      return 0;
    }
    else
      fprintf(stderr,"# Note: swapping byte-order!\n");
  }
  bytes += fread(&header->num_templates,1,4,fp);
  if (!byteorder_same) swap(4,1,(char *)&header->num_templates);
  bytes += fread(&header->f_low,1,4,fp);
  if (!byteorder_same) swap(4,1,(char *)&header->f_low);
  bytes += fread(&header->srate,1,4,fp);
  if (!byteorder_same) swap(4,1,(char *)&header->srate);
  header->mass = (struct MassPair *)realloc((void *)header->mass,
                                            8*(header->num_templates));
  if (header->mass==NULL) {
    GR_start_error("read_header()",rcsid,__FILE__,__LINE__);
    GR_report_error("Unable to allocate %d bytes of memory\n",
                    8*(header->num_templates));
    GR_end_error();
    return 0;
  }
  bytes += fread(header->mass,1,8*(header->num_templates),fp);
  if (!byteorder_same) swap(4,2*(header->num_templates),(char *)header->mass);
  bytes -= headerlength + 16 + 8*(header->num_templates);
  return !bytes;
}
#undef DEFAULT_HEADER_SIZE


/* Read signal files---see documentation */
int read_segment(FILE *fp, struct MF_Segment *segment, int num_templates)
{
  int bytes=0;
  bytes += fread(&segment->datastart,1,8,fp);
  if (!byteorder_same) swap(8,1,(char *)&segment->datastart);
  bytes += fread(&segment->is_gauss,1,4,fp);
  if (!byteorder_same) swap(4,1,(char *)&segment->is_gauss);
  segment->signal = (struct MF_Signal *)realloc((void *)segment->signal,
		                                32*num_templates);
  if (segment->signal==NULL) {
    GR_start_error("read_segment()",rcsid,__FILE__,__LINE__);
    GR_report_error("Unable to allocate %d bytes of memory\n",32*num_templates);
    GR_end_error();
    return 0;
  }
  bytes += fread(segment->signal,1,32*num_templates,fp);
  if (!byteorder_same) swap(4,8*num_templates,(char *)segment->signal);
  bytes -= 12 + 32*num_templates;
  return !bytes;
}


/* Usage message for the program */
void usage (char *program)
{
  fprintf(stderr,"Usage: %s [hl:]\n"
  "Options:\n"
  "  -h           print this message and exit\n"
  "  -l logfile   follows injection logfile\n",program);
}
