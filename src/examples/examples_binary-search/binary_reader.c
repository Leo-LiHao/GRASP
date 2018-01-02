/* GRASP: Copyright 1997,1998  Bruce Allen */
/* binary_reader.c
 * (based on read_multi.c by PRB, based on readmulti.c by JC)
 *
 * program to read the output from the binary_search program
 * running on 40m data
 *
 * this version by Jolien Creighton (e-mail: jolien@tapir.caltech.edu)
 * and Patrick Brady (e-mail: patrick@tapir.caltech.edu)
 * and Bruce Allen (e-mail ballen@dirac.phys.uwm.edu)
 * and R. Balasubramanian (e-mail bala@chandra.phys.uwm.edu)
 */


/* Examples:
 *
 * To print the SNR and segment number of each segment:
 *
 * binary_reader
 *
 * To print the segment start time and offset for every data segment that
 * passes an outlier test and has a maximum SNR greater than 10 where the
 * maximization is over all templates that have r^2 less than 2.5:
 *
 * binary_reader -o -t 10.0 -r 2.5 -f tO
 *
 * To print a useful help message
 *
 * binary_reader -h
 *
 */

#include "grasp.h"
#include <string.h>
#include <float.h>

extern char *optarg;
extern int optind;

const char *rcsid="$Id: binary_reader.c,v 1.14 1999/03/30 18:21:17 ballen Exp $";

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
  extern int getopt (int argc, char *const *argv, const char *shortopts);
  void usage(char *program);

  float stat_max=FLT_MAX,threshold=0;
  int labels=0,outlier_test=0,print_templates=0;
  int template_min=0,template_max=-1,median_maximization=0;
  int fmtlen=2,count=0;
  char fmt[16]="Ss",sep[4]="\t",*mfpath;
  struct MF_Header header={0};
  double timeofarrival,fsamp;
  int flag_compare=0;
  FILE *fp;
  int particularsegment=-1;
  char ps_file[40];
  int flag_ps=0;
  FILE *fp_ps;

  /* parse the optional arguments */
  while (1){
    int c;
      
    /* call standard C library option parser */
    c = getopt(argc, argv, "hmf:lor:t:L:T:Ms:c");

    if (c == -1 && optind == 1) { 
      /* If there are no arguments print to stdout */
      fprintf(stderr,"# SNR maximised over template bank and arrival time\n");
      break;
    } else if (c == -1)
      break;
      
    switch (c) {

    case 'h': /* Print a simple usage message */
      usage(argv[0]);
      exit(0);

    case 'm': /* Print a list of templates as:  No. m_1 m_2 */
      fprintf(stderr,"# Printing template list,  all other flags ignored\n");
      print_templates = 1;
      break;

    case 'f': /* Print format */
      fmtlen = strlen(optarg);
      if (fmtlen>16) {
	fprintf(stderr,"# Warning: print format truncated at 16 characters\n");
	fmtlen = 16;
      }
      strncpy(fmt,optarg,16);
      break;

    case 'l': /* Print field labels */
      labels = 1;
      sprintf(sep,"; ");
      break;

    case 'o': /* Reject data segments with outliers */
      outlier_test = 1;
      fprintf(stderr,"# Tested for outliers\n");
      break;

    case 'r': /* Exclude templates failing r^2 descriminant */
      stat_max = atof(optarg);
      fprintf(stderr,"# SNR if Bruce's r^2 descriminant > %f\n",stat_max);
      break;

    case 't': /* Exclude segments in which max snr less than threshold */
      threshold = atof(optarg);
      fprintf(stderr,"# SNR threshold: %f\n",threshold);
      break;

    case 'L': /* Exclude results from templates 0-->(optarg)-1 */
      template_min = atoi(optarg);
      fprintf(stderr,"# Excluded templates 0--%i\n",template_min-1);
      break;
	      
    case 'T': /* Print SNR from specific template---identified by number */
      template_min = template_max = atoi(optarg);
      fprintf(stderr,"# SNR for template number %i\n",template_max);
      break;

    case 'M': /* Maximize over median-renormalized SNR */
      median_maximization = 1;
      break;
        
    case 's': /* particular segment to be analysed */
        particularsegment = atoi(optarg);
        if(particularsegment<0){
            printf("segment number cannot be negative\n");
            exit(1);
        }
        flag_ps = 1;
        break;
    case 'c': /* to interpret the offsets as times of arrival */
 	flag_compare=1;
	break;
              
    default: /* If things fail for some reason */
      fprintf(stderr,"# getopt returned character code 0%o ??", c);
    }
  }
  
  /* Deal with possible errors */
  if (optind < argc) {
    fprintf(stderr,"# non-option ARGV-elements: ");
    while (optind < argc) fprintf(stderr,"%s ", argv[optind++]);
  }

  /* get path to read the data */
  mfpath = getenv("GRASP_MFPATH");
  if (mfpath==NULL) {
    GR_start_error("main",rcsid,__FILE__,__LINE__);
    GR_report_error("The environment variable GRASP_MFPATH must be set.\n");
    GR_end_error();
    return 1;
  }
  else
    fprintf(stderr,"# Using data from path GRASP_MFPATH = %s\n",mfpath);

  /* Read the header file */
  fp = fopen_segment(mfpath,-1);
  if (fp==NULL) {
    GR_start_error("main()",rcsid,__FILE__,__LINE__);
    GR_report_error("Unable to read header file: %s/signal.header\n",mfpath);
    GR_end_error();
    exit(1);
  }
  if (!read_header(fp,&header)) {
    GR_start_error("main()",rcsid,__FILE__,__LINE__);
    GR_report_error("Error reading header file: %s/signal.header\n",mfpath);
    GR_end_error();
    exit(1);
  }
  fclose(fp);
  fsamp = header.srate;

  if (template_max<0) template_max = header.num_templates - 1;

  /* If print_templates,  then print out the list of templates and return */
  if(print_templates){
      int i;
      fprintf(stdout,"Templates used to produce these signals. Format: No. m1 m2\n\n");
      for(i=template_min;i<=template_max;i++)
	  fprintf(stdout,"%i \t%f %f\n",i,header.mass[i].one,header.mass[i].two);
      return(0);
  }
  if(flag_ps==1){
      int template;
      if((fp=fopen_segment(mfpath,particularsegment))!=NULL){
          struct MF_Segment segment={0};
          sprintf(ps_file,"segment.%d",particularsegment);
          if((fp_ps = fopen(ps_file,"w"))==NULL){
              fprintf(stderr,"# Can't open file ./segment.%d\n",particularsegment);
              exit(1);
          }
          if (!read_segment(fp,&segment,header.num_templates)) {
              GR_start_error("main()",rcsid,__FILE__,__LINE__);
              GR_report_error("Unable to read segment file: %s/signal.%05d\n",mfpath,particularsegment);
              GR_end_error();
              exit(1);
          }          
          for (template=template_min;template<=template_max;template++){
	      if(flag_compare==1) 
                  fprintf(fp_ps,"%d %f %f %f %f %d\n",		
                          -segment.signal[template].offset,
                          segment.signal[template].snr_max,
                          segment.signal[template].snr_max_pow_renorm,
                          segment.signal[template].snr_max_med_renorm,
                          segment.signal[template].distance,
                          segment.signal[template].num_crossing);
              else
                  fprintf(fp_ps,"%d %f %f %f %f %d\n",		
                          segment.signal[template].offset,
                          segment.signal[template].snr_max,
                          segment.signal[template].snr_max_pow_renorm,
                          segment.signal[template].snr_max_med_renorm,
                          segment.signal[template].distance,
                          segment.signal[template].num_crossing);
          }
      }
      else{
          printf("Cannot Open the signal file for segment %d\n",particularsegment);
          exit(1);
      }
      fclose(fp_ps);
      fclose(fp);
  }
  while ((fp=fopen_segment(mfpath,count++))!=NULL) {

    struct MF_Segment segment={0};
    int max_template=0;
    int print_flag;

    /* Read data file */
    if (!read_segment(fp,&segment,header.num_templates)) {
      GR_start_error("main()",rcsid,__FILE__,__LINE__);
      GR_report_error("Unable to read segment %05d file: %s/signal.%05d\n",
		count-1,mfpath,count-1);
      GR_end_error();
      exit(1);
    }
    fclose(fp);

    print_flag = (!outlier_test || segment.is_gauss);

    if (print_flag) {
      int template;
      float max=threshold;
      for (template=template_min;template<=template_max;template++)
	if (segment.signal[template].stat<stat_max) {
	  float tmp;
          if (median_maximization)
	    tmp = segment.signal[template].snr_max_med_renorm;
	  else
            tmp = segment.signal[template].snr_max;
	  if (tmp>max) {
	    max_template = template;
            max = tmp;
	  }
	}
      print_flag &= (max > threshold);
    }

    if (print_flag) {
      int i;
      for (i=0;i<fmtlen;i++) {
	switch (fmt[i]) {
	case 't':
	  if (labels) printf("time ");
	  printf("%f",segment.datastart);
	  break;
	case 's':
	  if (labels) printf("segment ");
	  printf("%05d",count-1);
	  break;
	case 'o':
	  if (labels) printf("Outliers? ");
	  printf("%d",!segment.is_gauss);
	  break;
	case 'T':
	  if (labels) printf("max template = ");
	  printf("%d",max_template);
	  break;
	case 'M':
	  if (labels) printf("max template masses = ");
	  printf("%.3f",header.mass[max_template].one);
	  if (labels) printf(", ");
	  else printf("%s",sep);
	  printf("%.3f",header.mass[max_template].two);
	  break;
	case 'D':
	  if (labels) printf("max template dist = ");
	  printf("%.3f",segment.signal[max_template].distance);
	  break;
	case 'S':
	  if (labels) printf("max template snr = ");
	  printf("%.2f",segment.signal[max_template].snr_max);
	  break;
	case 'X':
	  if (labels) printf("max template snr (pow renorm) = ");
	  printf("%.2f",segment.signal[max_template].snr_max_pow_renorm);
	  break;
	case 'Y':
	  if (labels) printf("max template snr (med renorm) = ");
	  printf("%.2f",segment.signal[max_template].snr_max_med_renorm);
	  break;
	case 'O':
	  if (labels) printf("max template offset = ");
          if(flag_compare==1){
              if(segment.signal[max_template].offset>0){
                  printf("\nError  : -c option inconsistent with these segment files\n");
                  printf("         : Recorded offset is the impulse offset\n");
                  printf("         : Please run the program without the -c option\n");
                  exit(1);
              }
              timeofarrival = -1.*segment.signal[max_template].offset/fsamp +(double)segment.datastart;
              printf("%d (%f)",-segment.signal[max_template].offset,timeofarrival);
          }
          else{
              if(segment.signal[max_template].offset<0){
                  printf("\nError  : The -c option has to be used in the command line for printing the offset.\n");
                  printf("         : Recorded offset in these segment files corresponds to the\n");
                  printf("         : time of arrival.\n");
                  exit(1);
              }
              printf("%d",segment.signal[max_template].offset);
          }
	  break;
	case 'P':
	  if (labels) printf("max template phase = ");
	  printf("%.3f",segment.signal[max_template].phase);
	  break;
	case 'R':
	  if (labels) printf("max template stat = ");
	  printf("%.2f",segment.signal[max_template].stat);
	  break;
	case 'N':
	  if (labels) printf("max template ncross = ");
	  printf("%d",segment.signal[max_template].num_crossing);
	  break;
	default:
	  break;
	}
	printf("%s",sep);
      }
      printf("\n");
      fflush(stdout);
    }
    free(segment.signal);
  }

  /* Finished */
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
    		GR_report_error("It does not appear to be either big- or little-endian!\n");
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
  header->mass = (struct MassPair *)realloc((void *)header->mass,8*(header->num_templates));
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
  segment->signal = (struct MF_Signal *)realloc((void *)segment->signal,32*num_templates);
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
  fprintf(stderr,"Usage: %s [hmf:lor:t:]\n"
  "Options:\n"
  "  -h           print this message and exit\n"
  "  -s segno     print signals for all the templates for the segment number segno to file segment.segno\n"
  "  -c           interprets the offsets stored as times of arrival rather than impulse offsets\n"        
  "  -m           prints the list of templates to standard out\n"
  "  -f fmt       format fmt for fields to print\n"
  "  -l           print labels for fields\n"
  "  -o           only print segments without outliers\n"
  "  -r rmax      reject filters with r^2 > rmax in maximization\n"
  "  -t threshold only print segments with maximum SNR > threshold\n"
  "  -L lower     only maximize over templates numbers > lower\n"
  "  -T tempno    consider only template number tempno--do not maximize\n"
  "  -M           maximize over median-renormalized SNR\n"
  "\n"
  "Format: fmt = [tsoTMDSOPRN] (default fmt = Ss)\n"
  "  t:           print segment start time\n"
  "  s:           print segment number\n"
  "  o:           print 1 if segment has outliers, 0 otherwise\n"
  "  T:           print number of maximum template\n"
  "  M:           print masses of maximum template\n"
  "  D:           print distance of SNR=1 for maximum template\n"
  "  S:           print max SNR of maximum template\n"
  "  X:           print max SNR of maximum template with power renorm\n"
  "  Y:           print max SNR of maximum template with median renorm\n"
  "  O:           print max offset of maximum template\n"
  "  P:           print phase of maximum template\n"
  "  R:           print r^2 test value of maximum template\n"
  "  N:           print number of crossings of maximum template\n",program);
}
