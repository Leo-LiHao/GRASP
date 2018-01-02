#include "grasp.h"
#define SECPERDAY (3600*24)

static char *rcsid="$Id: utctime.c,v 1.6 1999/05/12 23:08:03 ballen Exp $\
\n$Name: RELEASE_1_9_8 $";

/* leap second table, from ftp://maia.usno.navy.mil/ser7/tai-utc.dat

 1972 JAN  1 =JD 2441317.5  TAI-UTC=  10.0
 1972 JUL  1 =JD 2441499.5  TAI-UTC=  11.0
 1973 JAN  1 =JD 2441683.5  TAI-UTC=  12.0
 1974 JAN  1 =JD 2442048.5  TAI-UTC=  13.0
 1975 JAN  1 =JD 2442413.5  TAI-UTC=  14.0
 1976 JAN  1 =JD 2442778.5  TAI-UTC=  15.0
 1977 JAN  1 =JD 2443144.5  TAI-UTC=  16.0
 1978 JAN  1 =JD 2443509.5  TAI-UTC=  17.0
 1979 JAN  1 =JD 2443874.5  TAI-UTC=  18.0
 1980 JAN  1 =JD 2444239.5  TAI-UTC=  19.0
 1981 JUL  1 =JD 2444786.5  TAI-UTC=  20.0
 1982 JUL  1 =JD 2445151.5  TAI-UTC=  21.0
 1983 JUL  1 =JD 2445516.5  TAI-UTC=  22.0
 1985 JUL  1 =JD 2446247.5  TAI-UTC=  23.0
 1988 JAN  1 =JD 2447161.5  TAI-UTC=  24.0
 1990 JAN  1 =JD 2447892.5  TAI-UTC=  25.0
 1991 JAN  1 =JD 2448257.5  TAI-UTC=  26.0
 1992 JUL  1 =JD 2448804.5  TAI-UTC=  27.0
 1993 JUL  1 =JD 2449169.5  TAI-UTC=  28.0
 1994 JUL  1 =JD 2449534.5  TAI-UTC=  29.0
 1996 JAN  1 =JD 2450083.5  TAI-UTC=  30.0
 1997 JUL  1 =JD 2450630.5  TAI-UTC=  31.0
 1999 JAN  1 =JD 2451179.5  TAI-UTC=  32.0
*/


/* Additional entries from
MJD = JD - 2400000.5

 1968 FEB  1 =JD 2439887.5  TAI-UTC=   4.2131700 S + (MJD - 39126.) X 0.002592 S
 1972 JAN  1 =JD 2441317.5  TAI-UTC=  10.0       S + (MJD - 41317.) X 0.0      S

 1970 JAN  1 =JD 2440587.5  TAI-UTC=  8.0
 1971 JAN 22 =JD 2440973.5  TAI-UTC=  9.0

*/


/* additional entries may be filled in following the same pattern, as */
/* future leap seconds are decided upon */

/*  see web site ftp://maia.usno.navy.mil/ser7/tai-utc.dat for updating table */
static time_t leaps[]={
		 ((2440587-2440587)*SECPERDAY),
		 ((2440973-2440587)*SECPERDAY),
		 ((2441317-2440587)*SECPERDAY),
		 ((2441499-2440587)*SECPERDAY),
                 ((2441683-2440587)*SECPERDAY),
		 ((2442048-2440587)*SECPERDAY),
		 ((2442413-2440587)*SECPERDAY),
		 ((2442778-2440587)*SECPERDAY),
		 ((2443144-2440587)*SECPERDAY),
		 ((2443509-2440587)*SECPERDAY),
		 ((2443874-2440587)*SECPERDAY),
		 ((2444239-2440587)*SECPERDAY),
		 ((2444786-2440587)*SECPERDAY),
		 ((2445151-2440587)*SECPERDAY),
		 ((2445516-2440587)*SECPERDAY),
		 ((2446247-2440587)*SECPERDAY),
		 ((2447161-2440587)*SECPERDAY),
		 ((2447892-2440587)*SECPERDAY),
		 ((2448257-2440587)*SECPERDAY),
		 ((2448804-2440587)*SECPERDAY),
		 ((2449169-2440587)*SECPERDAY),
		 ((2449534-2440587)*SECPERDAY),
		 ((2450083-2440587)*SECPERDAY),
		 ((2450630-2440587)*SECPERDAY),
		 ((2451179-2440587)*SECPERDAY)
		};

struct tm *utctime(const time_t *tp) {
	static struct tm returntm;
	struct tm *gmt;
	time_t inputtime;
	const int arraylen=sizeof(leaps)/sizeof(time_t);
	const int maxtested=946684823.0;
	int i=0;

	/* make sure that the argument is in-range */
	if (*tp<0) {
		GR_start_error("utctime()",rcsid,__FILE__,__LINE__);
		GR_report_error("The argument to utctime() must be a pointer to a time_t\n");
		GR_report_error("type in the range >=0. Input value was %d\n",(int)*tp);
		GR_report_error("This function does not have a table of Leap Sec before Jan 1, 1970.\n");
		GR_end_error();
	}

/* see web site ftp://maia.usno.navy.mil/ser7/tai-utc.dat for updating test boundary! */
	if (*tp>maxtested) {
		GR_start_error("utctime()",rcsid,__FILE__,__LINE__);
		GR_report_error("The argument to utctime() must be a pointer to a time_t\n");
		GR_report_error("type in the range < %d . Input value was %d\n",(int)maxtested,(int)*tp);
		GR_report_error("This function does not have a table of Leap Sec after Dec 31, 1999.\n");
		GR_report_error("It should be re-coded with this information, if possible.\n");
		GR_end_error();
	}

	/* This function assues that gmtime() is broken, check this */
	inputtime=(2441317-2440587)*SECPERDAY;
	gmt=gmtime(&inputtime);
	returntm=*gmt;
	/* since this is a leap year, if gmtime() NOT broken, sec field will be 60 */
	if (returntm.tm_sec==60) {
		GR_start_error("utctime()",rcsid,__FILE__,__LINE__);
		GR_report_error("Warning: On this machine, the standard C-library function gmtime() appears\n");
		GR_report_error("to correctly handle leap seconds.  We will use this instead of utctime()!\n");
		GR_end_error();
		return gmtime(tp);
	}

	/* make a local copy of the epoch time */
	inputtime=*tp;

	/* scan along table to see if it corresponds to any leap second since 1970 */
	while (i<arraylen && (leaps[i]+i-1)<inputtime)
		i++;

	/* if it does correspond to a leap second, then */
	if ((leaps[i]+i-1)==inputtime) {
		inputtime-=i;
		gmt=gmtime(&inputtime);
		returntm=*gmt;
		/* set the seconds field to 60, rather than 0,..,59 */
		returntm.tm_sec=60;
	} else {
		/* else skip the appropriate number of leap seconds */
		inputtime-=(i-1);
		gmt=gmtime(&inputtime);
		returntm=*gmt;
	}

	/* and return the UTC time structure */
	return &returntm;
}


struct tm *gpstime(const time_t *tp) {
	time_t inputtime;
	inputtime=*tp;

	/* shift the time to GPS time */
	inputtime+=UTCTOGPS;

	/* and return the appropriate structure */
	return utctime(&inputtime);
}

