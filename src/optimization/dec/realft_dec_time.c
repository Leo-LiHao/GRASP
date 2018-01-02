/* GRASP: Copyright 1997,1998  Bruce Allen */
/* This is a simple timing program to comnpare the speed of 
replacement realft routines. */

static char *rcsid="$Id: realft_dec_time.c,v 1.3 1998/01/23 17:59:40 ballen Exp $\n$Name: RELEASE_1_9_8 $";

#include <stdio.h>
#include <math.h>  
#include <string.h>

#define NREPEAT 25
#define NPOINT1 524288
#define NPOINT2 262144 

void realft(float *,unsigned long,int); 
void printvec(char *title,float *f, int n);

main()
{
	float x,norm,*array1,*array2,*work1,*work2;
	int i;


	array1=(float *)malloc(sizeof(float)*NPOINT1);
	work1=(float *)malloc(sizeof(float)*NPOINT1);

	norm=2.0/NPOINT1;

	for (i=0;i<NPOINT1;i++) {
          x=0.0000001*i;
	  array1[i]=cos(x)+x*x-6.0*x+5.0;}


	for (i=0;i<NREPEAT;i++){
	memcpy(work1,array1,sizeof(float)*NPOINT1);  
	realft(work1-1,NPOINT1,1); 
 	realft(work1-1,NPOINT1,-1);
        } 


        array2=(float *)malloc(sizeof(float)*NPOINT2);
	work2=(float *)malloc(sizeof(float)*NPOINT2);
	norm=2.0/NPOINT2;
	
        for (i=0;i<NPOINT2;i++) {
          x=0.0000001*i;
	  array2[i]=cos(x)+x*x-6.0*x+5.0;}

	
	for (i=0;i<NREPEAT;i++){
	memcpy(work2,array2,sizeof(float)*NPOINT2);  
	realft(work2-1,NPOINT2,1); 
 	realft(work2-1,NPOINT2,-1); 
	} 

	return 0;
}
