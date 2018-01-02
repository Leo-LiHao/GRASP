/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"

int main() {
	struct Scope Grid;

	/* Set parameters for the inspiral search in CIT 40 meter */
	Grid.m_mn=0.9;
	Grid.m_mx=3.0;
	Grid.theta=0.307;
	Grid.dp=2*0.008;
	Grid.dq=2*0.0006;
	Grid.f_start=140.0;

	/* construct template set covering parameter space */
	template_grid(&Grid);

	/* create a postscript file showing locations of templates */
	plot_template("templates_40meter.ps",Grid,15,1);
	return 0;
}
