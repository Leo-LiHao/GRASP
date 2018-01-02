/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"

int main(int argc, char **argv)
{
  struct cubic_grid grid;

  /* Set grid parameters. */
  grid.n=13;
  grid.m_mn=0.8;
  grid.m_mx=3.2;
  grid.match=0.98;
  grid.angle=0.0;
  grid.order=4;
  grid.srate=50000.0;
  grid.flo=120.0;
  grid.ftau=140.0;
  grid.detector=15; /* Smooth fit to Caltech 40m prototype. */
  /* grid.detector=8;  Caltech 40m prototype. */
  /* grid.detector=1;  LIGO initial interferometer. */
  /* grid.detector=12; LIGO advanced interferometer. */

  /* Generate grid of cubic-fit coefficients */
  generate_cubic(grid,"detectors.dat",
		 "cubic_coef_40meter_m=0.8-3.2.ascii",
		 "cubic_coef_40meter_m=0.8-3.2.log");
  return 0;
}
