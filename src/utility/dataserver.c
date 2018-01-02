/* GRASP: Copyright 1997,1998  Bruce Allen */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <strings.h>

static char *rcsid="$Id: dataserver.c,v 1.3 1998/01/23 17:59:52 ballen Exp $\
\n$Name: RELEASE_1_9_8 $";

int      PF_urlopen_timeout  (char*, char*, int*, int);
int      PF_urlopen          (char*, char*, int*);
void       PF_urlopen_perror   (int);

extern char *line;
extern int info;

#define SILENT  0
#define QUIET   1
#define NORMAL  2
#define CHATTY  3
#define VERBOSE 4
#define DEBUG   5


#define ERRNO_NOT_SET                  0

/* normal return */
#define HTTP_CONTENT_READ              1

#define URL_IS_NOT_HTTP                2
#define URL_IS_NULL_STRING             3
#define TOO_MANY_HEADER_LINES          4
#define TIMED_OUT                      5
#define NO_DNS_ENTRY                   6
#define BIND_FAILED                    7
#define CONNECT_FAILED                 8
#define WRITE_GET_REQUEST_FAILED       9

/* add this to error-code from URL header */
#define HTTP_ERROR                     100

#include <stdio.h>
int urlopen_timeout(char *, int);
void urlopen_perror(FILE *);



/***************************************************************************/
/*    Program to illustrate extraction of data from CACR HPSS using        */
/*    metadata and data web servers.                                       */
/*    Roy Williams June 1997                                               */
/***************************************************************************/


FILE *fp_metadata;

/* First call this function saying the start and end times for the
   data files that you want. The times are Unix time, ie seconds 
   since 1/1/70
*/

int init_fd_iterator(long Tstart, long Tend)
{
	char s[1024], mime_type[1024];
	int size, fd_metadata;

	GR_start_error("init_fd_iterator()",rcsid,__FILE__,__LINE__);
	GR_report_error("Dumping frames from Tstart = %ld (%s)\n", Tstart, ctime(&Tstart));
	GR_report_error("                 to Tend   = %ld (%s)\n", Tend, ctime(&Tend));
	GR_end_error();

	sprintf(s,
	"http://sara.cacr.caltech.edu/ligo/getframes.htf?Tstart=%ld&Tend=%ld",
	Tstart, Tend);

	GR_start_error("init_fd_iterator()",rcsid,__FILE__,__LINE__);
	GR_report_error("Opening %s\n", s);
	GR_end_error();

	fd_metadata = PF_urlopen_timeout(s, mime_type, &size, 30);
	if(fd_metadata < 0){
		PF_urlopen_perror(0);
		return -1;
	}
	fp_metadata = fdopen(fd_metadata, "r");
	return 0;
}

/* Now repeatedly call this function to get a fil descriptor for the next
   frame file, until it returns -1, signifying the end
*/
int filedes()
{
	char line[1024], tag[1024], file[1024], nframestr[1024];
	char s[1024], mime_type[1024];
	int size, fd_data;
	static int first=1;
	char *tstring;
	long tstart,tend;
	

	if (first) {
		first=0;
		tstring=getenv("GRASP_FRAMETS");
		if (tstring==NULL) {
			GR_start_error("filedes()",rcsid,__FILE__,__LINE__);
			GR_report_error("GRASP: filedes(): the environment variable GRASP_FRAMETS is not defined!\n");
			GR_report_error("GRASP: This should be the integer number of seconds after Jan 1 1970 for first frame.\n");
			GR_end_error();
			abort();
		}
		tstart=atol(tstring);

		tstring=getenv("GRASP_FRAMETE");
		if (tstring==NULL) {
			GR_start_error("filedes()",rcsid,__FILE__,__LINE__);
			GR_report_error("GRASP: filedes(): the environment variable GRASP_FRAMETE is not defined!\n");
			GR_report_error("GRASP: This should be the integer number of seconds after Jan 1 1970 for last frame.\n");
			GR_end_error();
			abort();
		}
		tend=atol(tstring);
		
		init_fd_iterator(tstart, tend);
	}

	if(fp_metadata == NULL) return -1;

	while(fgets(line, 256, fp_metadata)){
		sscanf(line, "%s %s %s", tag, file, nframestr);
		if(strcmp(tag, "File:") != 0) continue;

		GR_start_error("filedes()",rcsid,__FILE__,__LINE__);
		GR_report_error("File is %s, Nframe is %s\n", file, nframestr);
		GR_end_error();

		strcpy(s, "http://bitty.cacr.caltech.edu:1080/cgi-bin/ligo_frame_server?");
		strcat(s, "FORMAT=raw&FILE=");
		strcat(s, file);

/*   Lets try this as SDF some time.... */
/*		strcat(s, "FORMAT=SDF&FILE="); */

		fd_data = PF_urlopen(s, mime_type, &size);
		if(fd_data < 0){
			PF_urlopen_perror(0);
			return -1;
		} else {
			GR_start_error("filedes()",rcsid,__FILE__,__LINE__);
			GR_report_error("Mime Type %s\n", mime_type);
			GR_end_error();
			return fd_data;
		}
	}
	fclose(fp_metadata);
	return -1;
}

/* Now repeatedly call this function to get a fil descriptor for the next
   frame file, until it returns -1, signifying the end
*/
int fd_iterator()
{
	char line[1024], tag[1024], file[1024], nframestr[1024];
	char s[1024], mime_type[1024];
	int size, fd_data;

	if(fp_metadata == NULL) return -1;

	while(fgets(line, 256, fp_metadata)){
		sscanf(line, "%s %s %s", tag, file, nframestr);
		if(strcmp(tag, "File:") != 0) continue;

		GR_start_error("fd_iterator",rcsid,__FILE__,__LINE__);
		GR_report_error("File is %s, Nframe is %s\n", file, nframestr);
		GR_end_error();

		strcpy(s, "http://bitty.cacr.caltech.edu:1080/cgi-bin/ligo_frame_server?");
		strcat(s, "FORMAT=raw&FILE=");
		strcat(s, file);

/*   Lets try this as SDF some time.... */
/*		strcat(s, "FORMAT=SDF&FILE="); */

		fd_data = PF_urlopen(s, mime_type, &size);
		if(fd_data < 0){
			PF_urlopen_perror(0);
			return -1;
		} else {
			GR_start_error("fd_iterator",rcsid,__FILE__,__LINE__);
			GR_report_error("Mime Type %s\n", mime_type);
			GR_end_error();
			return fd_data;
		}
	}
	fclose(fp_metadata);
	return -1;
}

#ifdef STANDALONE 

int main(int argc, char** argv)
{
	char *buff;
	long buffSize;
	long Tstart, Tend;
	int fd_data, framenum;
	char file[1024];
	struct FrFile *iFile;
	struct FrameH *frame;

/* Set up the frame library */
	buffSize = 100000;
	buff = malloc(buffSize);
	if(FrLibIni(NULL,stderr,0) == NULL) {
		GR_start_error("main() example",rcsid,__FILE__,__LINE__);
		GR_report_error("Error during initialisation\n"
		    "  Last errors are:\n%s",FrErrorGetHistory());
		GR_end_error();
		exit(1);
	}

/* Figure out the times the user wants data */
/* Defaults are Tstart = Tue Nov 15 06:13:20 1994
				Tend   = Mon Nov 21 03:24:08 1994
*/
	if(argc > 1) Tstart = atoi(argv[1]); else Tstart = 784880000;
	if(argc > 2) Tend   = atoi(argv[2]); else Tend =   785388248;

/* Start the fd iterator */
	init_fd_iterator(Tstart, Tend);

/* Loop over the iterator */
	while((fd_data = fd_iterator()) >= 0){

/* Say something intelligent about the frames in this file */
		iFile = FrFileINewFd(fd_data,buff,buffSize);
		if(iFile == NULL) {
			GR_start_error("main() example",rcsid,__FILE__,__LINE__);
			GR_report_error("Error during file opening: Last errors are:\n%s",
				FrErrorGetHistory());
			GR_end_error();
			exit(1);
		}

		GR_start_error("main() example",rcsid,__FILE__,__LINE__);
		GR_report_error("Reading Frames...\n");
		GR_end_error();
		framenum = 0;
		while((frame = FrameRead(iFile)) != NULL) {
			GR_start_error("main() example",rcsid,__FILE__,__LINE__);
			GR_report_error("%s %2d %3d %6ld %10ld %10ld %10.6lf\n",
				file, 
				++framenum,
				frame->run, 
				frame->frame,
				frame->UTimeS, 
				frame->UTimeN, 
				frame->dt);
			GR_end_error();
			FrameFree(frame);
		}
		GR_start_error("main() example",rcsid,__FILE__,__LINE__);
		GR_report_error("Closing File...\n");
		GR_end_error();

		FrFileOEnd(iFile);
	}
	return(0);
}
#endif

/* This stuff is hack it is in the Paraflow library */
#include <stdarg.h>
void PF_print(int infolevel, char* fmt, ...)
{
    char line[256];
    va_list alist;
    va_start(alist, fmt);
    vsprintf(line, fmt, alist);
    va_end(alist);
    GR_start_error("PF_print()",rcsid,__FILE__,__LINE__);
    GR_report_error("%s", line);
    GR_end_error();
 
}

void _PF_out(int info, char* out)
{
    GR_start_error("PF_out()",rcsid,__FILE__,__LINE__);
    GR_report_error("%s", out);
    GR_end_error();
}


/************************************************************************/
/*                              urlopen.c                               */
/*       To turn a URL on a remote machine into a file pointer          */
/*            Mark Niedengard, Roy Williams, George Kremenek            */
/************************************************************************/

/* must link with -lnsl -lucb -lsocket on Solaris */

#include <setjmp.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <stdio.h>
#include <netinet/in.h>
#include <netdb.h>
#include <unistd.h>
#include <signal.h>


/***************************TEST TEST *********************************/
#ifdef STANDALONE1
main(int argc, char **argv)
/* simple test program */
{
	char sTmp[512], mime_type[256];
	int c, size;
	FILE *fp;
	int fd;

    strcpy(sTmp, argv[1]);

	if((fd = PF_urlopen_timeout(sTmp, mime_type, &size, 30)) >= 0){
		GR_start_error("main test",rcsid,__FILE__,__LINE__);
		GR_report_error("Content-type is %s\nContent-length is %d\n", mime_type, size);
		GR_end_error();
		fp = fdopen(fd, "r");
		while((c = getc(fp)) != EOF) putchar(c);
	} else {
		PF_urlopen_perror(NORMAL);
	}
}
#endif
/***************************TEST TEST *********************************/

int debug = 0;
int PF_urlopen_errno = ERRNO_NOT_SET;

/* globals */
jmp_buf   my_environment;

void
catch(int signo)
/* Signal catcher for the alarm so we can error on timeout */
{
   signal(SIGALRM, catch);
  /* caught allarm */
  /* abort I/O operation in progres by non-local jump (dirty trick!) */
  longjmp(my_environment, 1);
}

int getURLbyParts(char* machine, int port, char* file )
/* Get the url once it is split into pieces */
{
	struct hostent *he;
	struct servent *se;
	struct sockaddr_in sin;
	int sock;
	char *getString;

	if(debug) {
		GR_start_error("getURLbyParts()",rcsid,__FILE__,__LINE__);
		GR_report_error("Machine %s, port %d, file %s\n", machine, port, file);
		GR_end_error();
	}

	getString = (char *)malloc(30 + strlen(file));
	assert(getString);

	he = gethostbyname( machine );
	if(!he){
		PF_urlopen_errno = NO_DNS_ENTRY;
		return -1;
	}

	sock = socket(he->h_addrtype, SOCK_STREAM, 0);

	bzero((caddr_t)&sin, sizeof(sin));
	sin.sin_family = he->h_addrtype;
	if (bind(sock, &sin, sizeof(sin)) < 0){
		PF_urlopen_errno = BIND_FAILED;
		return -1;
	}
	bcopy(he->h_addr, &sin.sin_addr, he->h_length);
	se = getservbyname("telnet", "tcp");

	sin.sin_port = htons(port);
	if (connect(sock, &sin, sizeof(sin)) < 0){
		PF_urlopen_errno = CONNECT_FAILED;
		return -1;
	}

	/* send GET message to http */
	sprintf(getString,"GET %s HTTP/1.0\r\n\r\n", file);
	if ( write( sock, getString, strlen(getString)) != strlen(getString)){
		PF_urlopen_errno = WRITE_GET_REQUEST_FAILED;
		return -1;
	}

	free (getString);
	return sock;
}

int getURL(char* url )
/* split the url into pieces then call getURLbyParts */
{
	char* s;
	char* machine;
	int machineLen; /* length of machine name */
	int port;
	char* file = 0;
	int fd;

	char http[] = "http://";
	int httpLen = strlen( http );
	if ( !url ) {
		PF_urlopen_errno = URL_IS_NULL_STRING;
		return -1;
	}
	if ( strncmp( http, url, httpLen) ) {
		PF_urlopen_errno = URL_IS_NOT_HTTP;
		return -1;
	}

	/* get the machine name */
	for( s = url + strlen(http) ; 
		*s && (*s != ':') && (*s != '/'); s++ );
	machineLen = s - url;
	machineLen -= httpLen;
	machine = (char*)malloc( machineLen + 1 ); 
	assert( machine );
	strncpy( machine, url + httpLen, machineLen );

/* this line added from original Niedengard version */
	machine[machineLen] = '\0';

	/* get port number */
	if ( *s == ':' ) {
		port = atoi( ++s );
		while( *s && (*s != '/') )
			s++;
	} else {
		port = 80;
	}

	/* get the file name */
	if ( *s == '/' )
		file = s;
	else
		file = "/";

	fd = getURLbyParts(machine, port, file );
	free( machine );
	return fd;
}

char *fgets_fd(char* buf, int nmax, int fd)
{
    int n = 0;
    int nbyte;
	if(fd < 0) 
		return NULL;
    while(1) {
        nbyte = read(fd, &buf[n], 1);
        if(nbyte == 0) break;
        n++;
        if(buf[n-1] == '\n' || n+1 >= nmax) break;
    }
    buf[n] = '\0';
    if(debug) {
	GR_start_error("fgets_fd()",rcsid,__FILE__,__LINE__);
	GR_report_error("strlen = %d, line is %s\n", n, buf);
	GR_end_error();

    }
    return buf;
}

#define SEE_HEADER

int PF_urlopen(char *url, char* mime_type, int *size)
/* This reads in the data from the socket and tosses out the header
   maybe we should pay more attention to the header.... */
{
	int fd;
	char *ptr;
	char line[256], out[256];
	int error_code;
	int headerlines = 0;
	int i;
	*size = 0;

	PF_print(DEBUG, "urlopen(%s) starts\n", url);

	fd = getURL(url );
	if(fd < 0) {
		return -1;
	}
	PF_print(DEBUG, "Reading from URL\n");
#define NL 10
#define CR 13
	fgets_fd(line, 256, fd);
	ptr = line;
	while(*ptr++ != ' ');
	sscanf(ptr, "%d", &error_code);
	while(*ptr++ != ' ');
	if(error_code != 200){
		PF_urlopen_errno = HTTP_ERROR + error_code;
		return -1;
	} else {
		while(fgets_fd(line, 256, fd)){
			if(line[0] == NL || (line[0] == CR && line[1] == NL)) break;
#ifdef SEE_HEADER
			if(headerlines++ > 8){
				PF_urlopen_errno = TOO_MANY_HEADER_LINES;
				return -1;
			}
			ptr = line; 
			while(*ptr != ':') ptr++; 
			*ptr++ = '\0';

			sprintf(out, "Header: %s then %s", line, ptr);
			for(i=0; i<strlen(line); i++)
				sprintf(out, "%2d ", line[i]); 
			sprintf(out, "\n");
			_PF_out(DEBUG, out);

			if(strcmp("Content-length", line) == 0
			|| strcmp("Content-Length", line) == 0)
				sscanf(ptr, "%d", size);
			if(strcmp("Content-type", line) == 0
			|| strcmp("Content-Type", line) == 0)
				sscanf(ptr, "%s", mime_type);
#endif
		}
		PF_urlopen_errno = HTTP_CONTENT_READ;
		return fd;
	}
}

int PF_urlopen_timeout(char* sTmp, char* mime_type, int* size, int wait) 
/* open the url with a timeout (seconds) returning error */
{
	int fd;

	if (setjmp(my_environment) == 0)  {

	  /* get here on direct call */
	  /* do not modify register vars in this leg of code */

	  signal(SIGALRM, catch);
	  alarm(wait); /* wait time in seconds before giving up */

	  fd = PF_urlopen(sTmp, mime_type, size);

      alarm(0); /* cancel alarm */

	  return(fd); /* OK */

	} else {
	  /* get here by calling longjmp */
	  PF_urlopen_errno = TIMED_OUT;
	  return -1;
	}
}

void PF_urlopen_perror(int infolevel)
/* Provides the reason for the prvious open attempt failing */
{
	switch(PF_urlopen_errno){
	case ERRNO_NOT_SET:
		PF_print(infolevel,"ERRNO_NOT_SET\n"); break;
	case URL_IS_NOT_HTTP :
		PF_print(infolevel,"URL_IS_NOT_HTTP\n"); break;
	case URL_IS_NULL_STRING :
		PF_print(infolevel,"URL_IS_NULL_STRING\n"); break;
	case TOO_MANY_HEADER_LINES :
		PF_print(infolevel,"TOO_MANY_HEADER_LINES\n"); break;
	case TIMED_OUT :
		PF_print(infolevel,"TIMED_OUT\n"); break;
	case NO_DNS_ENTRY :
		PF_print(infolevel,"NO_DNS_ENTRY\n"); break;
	case BIND_FAILED :
		PF_print(infolevel,"BIND_FAILED\n"); break;
	case CONNECT_FAILED :
		PF_print(infolevel,"CONNECT_FAILED\n"); break;
	case WRITE_GET_REQUEST_FAILED :
		PF_print(infolevel,"WRITE_GET_REQUEST_FAILED\n"); break;
	}

	if(PF_urlopen_errno >= HTTP_ERROR)
		PF_print(infolevel,"HTTP_ERROR %d\n", PF_urlopen_errno-HTTP_ERROR);
}



#include <stdio.h>
#include <fcntl.h>

#ifdef STANDALONE
unsigned char buf[4096];

main()
{
        int fd;

        fd = open("junk", O_RDONLY);
        bread(fd, buf, 16);
        bread(fd, buf, 1);
        bread(fd, buf, 170);
        bread(fd, buf, 1000);
        bread(fd, buf, 1);
        bread(fd, buf, 16);
        bread(fd, buf, 6);
        close(fd);
}
#endif

bread(int fd, char* buf, int nbytew)
{
        int i, got, nbyte;
        nbyte = 0;
        do {
                got = _bread(fd, buf+nbyte, nbytew-nbyte);
                nbyte += got;
        } while(got > 0 && nbyte < nbytew);
        return nbyte;
}

int    first      = 1;
char*  bread_buf;
int    maxnbuf    = 32*1024*1024;
int    nbuf;
int    ibuf;
int    bread_debug= 0;

breset()
{
        ibuf = 0;
}

_bread(int fd, char* requestbuf, int requestsize)
{
        int ncopy;
        if(first){
                while((bread_buf = (char*)malloc(maxnbuf)) == NULL)
                        maxnbuf /= 2;
                if(bread_debug) printf("bread: %d bytes in buffer\n",
maxnbuf);
                ibuf = nbuf = 0;
                first = 0;
        }
        if(ibuf == nbuf){
                nbuf = read(fd, bread_buf, maxnbuf);
                if(bread_debug) printf("bread: got %d bytes (asked for
%d)\n", 
                        nbuf, maxnbuf);
                if(nbuf == 0) return 0;
                ibuf = 0;
        }
        ncopy = nbuf - ibuf;
        if(requestsize < ncopy) ncopy = requestsize;
        if(bread_debug > 1) 
                printf("bread: copying %d (of %d) from %d\n", ncopy,
requestsize, ibuf);
        bcopy(bread_buf + ibuf, requestbuf, ncopy);
        ibuf += ncopy;
        return ncopy;
}

