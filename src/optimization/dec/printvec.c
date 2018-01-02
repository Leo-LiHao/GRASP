/* GRASP: Copyright 1997,1998  Bruce Allen */
#include <stdio.h>

void printvec(char *title,float f[], int n)

{
	int i;
	printf("%s\n",title);
	for(i=0; i<n; i++)
		printf("%d\t%f\n",i,f[i]);
	printf("\n");
}

