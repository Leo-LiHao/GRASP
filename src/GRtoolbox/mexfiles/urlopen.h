
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

