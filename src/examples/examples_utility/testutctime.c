#include <time.h>
#include <stdio.h>
#include <grasp.h>
struct tm *utctime(const time_t *tp);
struct tm *gpstime(const time_t *tp);

void printone(time_t time1) {
        time_t gpstimet;
        char utc[64],gmt[64],gps[64];

	/* calculate GPS time */
        gpstimet=time1-UTCTOGPS;

	/* construct character strings with appropriate time stamps */
        strcpy(gmt,asctime(gmtime(&time1)));
        strcpy(utc,asctime(utctime(&time1)));
        strcpy(gps,asctime(gpstime(&gpstimet)));

	/* print... */
	printf("C-time: %10d GPS Time %10d  UTC: %.24s GPS: %.24s gmtime: %s",(int)time1,(int)gpstimet,utc,gps,gmt);
	return;
}

void testutctime(time_t time1) {
	time1--;
        printone(time1++);
        printone(time1++);
        printone(time1++);
        printf("\n");
	return;
}

int main() {
	testutctime(1);
	testutctime(33350400);
	testutctime(63072000);
	testutctime(78796801);
	testutctime(94694402);
	testutctime(126230403);
	testutctime(315964811);
	testutctime(784880277);
	testutctime(911110677);
	return 0;
}
