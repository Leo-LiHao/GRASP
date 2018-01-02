/* Compile with mex -L/usr/lib -lc mxUrlopen.c                                  */
/*										*/
/* Steve Drasco 								*/
/* Summer 1998									*/

#include "mex.h"
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
#include <stdarg.h>
#include "urlopen.h"
#include "info.h"
int      PF_urlopen_timeout  (char*, char*, int*, int);
int      PF_urlopen          (char*, char*, int*);
void     PF_urlopen_perror   (int);

void urlgetput(char *url, char *path)
{
	FILE	*fp, *fpextern;
	char	mime_type[256];
        int	c, size, fd;
	
	fp = fopen(path,"w");	
	
	/* hacked this out of Roy's test function on top of urlopen.c */

        if((fd = PF_urlopen(url, mime_type, &size)) >= 0){
                fprintf(stderr, "Content-type is %s\nContent-length is %d\n", mime_type, size);
                fpextern = fdopen(fd, "r");
                /* while((c = getc(fpextern)) != EOF) putc(c, fp); */
        } else {
                PF_urlopen_perror(NORMAL);
        }

	/* end of hacking */

	fclose(fp);

}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	char	*url, *path;
	int	status, urlLength, pathLength;
	
        /* mex checks */
	if (nrhs != 2)
                mexErrMsgTxt("Need 2 inputs ('url','filepath/filename')");
	if (nlhs != 0)
		mexErrMsgTxt("No out put arguments");
	if (mxIsChar(prhs[0]) != 1 || mxIsChar(prhs[1]) != 1)
		mexErrMsgTxt("Both inputs must be strings");

	/* allocate memory for strings */
	urlLength = mxGetN(prhs[0]) + 1;
	pathLength = mxGetN(prhs[1]) + 1;
	url = mxCalloc(urlLength, sizeof(char));
	path = mxCalloc(pathLength, sizeof(char));

	/* copy matlab strings to C strings */
	status = mxGetString(prhs[0], url, urlLength);
	if( status != 0)
		mexWarnMsgTxt("Not enough space. URL string is truncated.");
	status = mxGetString(prhs[1], path, pathLength);
        if( status != 0)
                mexWarnMsgTxt("Not enough space. File path string is truncated.");

	/* call url function here */
		urlgetput(url, path);	
}
 

/************************************************************************/
/*                              urlopen.c                               */
/*       To turn a URL on a remote machine into a file pointer          */
/*            Mark Niedengard, Roy Williams, George Kremenek            */
/************************************************************************/

void PF_print(int infolevel, char* fmt, ...)
{
        char line[256];
    va_list alist;
        va_start(alist, fmt);
        vsprintf(line, fmt, alist);
        va_end(alist);
        fprintf(stderr, "%s", line);
}
void _PF_out(int info, char* out)
{
        fprintf(stderr, "%s", out);
}

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

        if(debug) printf("Machine %s, port %d, file %s\n", machine, port, file);

        getString = (char *)mxMalloc(30 + strlen(file));
        assert(getString);

        he = gethostbyname( machine );
        if(!he){
                PF_urlopen_errno = NO_DNS_ENTRY;
                return -1;
        }

        sock = socket(he->h_addrtype, SOCK_STREAM, 0);

        bzero((caddr_t)&sin, sizeof(sin));
        sin.sin_family = he->h_addrtype;
/*      Modified Bruce Allen 7/6/99 */
/*      if (bind(sock, &sin, sizeof(sin)) < 0){ */
        if (bind(sock, (struct sockaddr *)(&sin), sizeof(sin)) < 0){
                PF_urlopen_errno = BIND_FAILED;
                return -1;
        }
        bcopy(he->h_addr, &sin.sin_addr, he->h_length);
        se = getservbyname("telnet", "tcp");

        sin.sin_port = htons(port);
/*      Modified Bruce Allen 7/6/99 */
/*      if (connect(sock, &sin, sizeof(sin)) < 0){ */
        if (connect(sock, (struct sockaddr *)(&sin), sizeof(sin)) < 0){
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
        machine = (char*)mxMalloc( machineLen + 1 );
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
    if(debug) printf("strlen = %d, line is %s\n", n, buf);
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

