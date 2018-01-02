/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"

#define LONGITUDE 118.133  /* degrees West: Caltech */

int main(int argc, char *argv[])
{
  time_t now;
  double lst_h,lst_m,lst_s;
  double gst_h,gst_m,gst_s;
  float lst,gst,longitude=LONGITUDE;
  char loc[128],utc[128];

  if (argc==2) now = (time_t)atoi(argv[1]);
  else time(&now);

  strftime(loc,128,"%H:%M:%S %Z %a %d %b %Y",localtime(&now));
  strftime(utc,128,"%H:%M:%S UTC %a %d %b %Y",gmtime(&now));
  lst = local_sidereal_time(now,longitude);
  gst = local_sidereal_time(now,0);
  lst_s = 60*modf(60*modf(lst,&lst_h),&lst_m);
  gst_s = 60*modf(60*modf(gst,&gst_h),&gst_m);

  printf("%ld seconds since 0h UTC 1 Jan 1970\n%s\n%s\n",(long)now,loc,utc);
  printf("%02u:%02u:%02u (%09.6f h) Local Sidereal Time\n",
	 (char)lst_h,(char)lst_m,(char)lst_s,lst);
  printf("%02u:%02u:%02u (%09.6f h) Greenwitch Sidereal Time\n",
	 (char)gst_h,(char)gst_m,(char)gst_s,gst);

  return 0;
}
