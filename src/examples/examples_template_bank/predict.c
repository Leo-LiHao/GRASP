/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"

void compare(char *name,struct Scope *Grid) {
	float area,tarea,ratio;
	int predict;

	/* construct template set covering parameter space */
	template_grid(Grid);

	/* calculate area of parameter space, divide by area of template */
	area=template_area(Grid);
	tarea=0.25*M_PI*(Grid->dp)*(Grid->dq);
	ratio=Grid->n_tmplt*tarea/area;
	predict=(int)(area/tarea);

	/* see how closely this reflects actual number of templates needed */
	printf("Case: %s needed %d templates (predicted %d, ratio %f)\n",
		name,Grid->n_tmplt,predict,ratio);
	return;
}

int main() {
	struct Scope Grid;
	char what[64];

	strcpy(what,"CIT 40-meter (Oct), Sathyaprakash, M=0.8->50");
	Grid.m_mn=0.8;
	Grid.m_mx=50.0;
	Grid.theta=0.978;
	Grid.dp=2.0*0.00227;
	Grid.dq=2.0*0.0352;
	Grid.f_start=120.0;
	compare(what,&Grid);

	strcpy(what,"CIT 40-meter (Nov - 120 Hz cut), Sathyaprakash, M=0.5->50");
	Grid.m_mn=0.5;
	Grid.m_mx=50.0;
	Grid.theta=1.025;
	Grid.dp=2.0*0.00255;
	Grid.dq=2.0*0.0332;
	Grid.f_start=120.0;
	compare(what,&Grid);

	strcpy(what,"CIT 40-meter (Nov - 140 Hz cut), Sathyaprakash, M=0.5->50");
	Grid.m_mn=0.5;
	Grid.m_mx=50.0;
	Grid.theta=0.964;
	Grid.dp=2.0*0.00213;
	Grid.dq=2.0*0.0320;
	Grid.f_start=140.0;
	compare(what,&Grid);

	strcpy(what,"Initial LIGO, Owen, M=1.0->500");
	Grid.m_mn=1.0;
	Grid.m_mx=500.0;
	Grid.theta=0.5066;
	Grid.dp=2.0*0.000162324;
	Grid.dq=2.0*0.00210929;
	Grid.f_start=140.0;
	compare(what,&Grid);
	plot_template("temp_list.ps",Grid,30,1);

	strcpy(what,"Initial LIGO, Owen, M=0.2->500");
	Grid.m_mn=0.2;
	compare(what,&Grid);

	strcpy(what,"Advanced LIGO, Owen, M=1.0->500");
	Grid.m_mn=1.0;
	compare(what,&Grid);
/*
	strcpy(what,"Advanced LIGO, Owen, M=0.2->500");
	Grid.m_mn=0.2;
	Grid.theta=0.4524;
	Grid.dp=2.0*0.000352231;
	Grid.dq=2.0*0.00396995;
	Grid.f_start=70.0;
	compare(what,&Grid);

*/

	return 0;
}
