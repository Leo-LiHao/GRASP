/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"

int main() {
	struct Scope Grid;
	float area;

	/* Specify the parameter space */
	Grid.m_mn=0.8;
	Grid.m_mx=50.0;
	Grid.f_start=140.0;

	/* find area of parameter space */
	area=template_area(&Grid);

	/* and print it */
	printf("The area in parameter space is %f seconds^2.\n",area);
	return 0;
}
