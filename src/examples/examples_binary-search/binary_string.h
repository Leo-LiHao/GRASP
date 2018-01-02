/* GRASP: Copyright 1997,1998  Bruce Allen */
/* Description of this data run */
const char *description=
"Program binary_search created on " __DATE__ " " __TIME__ "\n"
"with parameters:\n\n"
#include "binary_params.ascii"
"\n";

/* Comment describing the signal file format */
const char *comment=
"Signal data output format:\n\n"
"The data format for the optimal filtering output is as follows:\n\n"
"W bytes of ascii text - a descriptor of what follows in free format.\n"
"W is encoded in the first line of this file (if not present, W==4096).\n"
"4 bytes integer 1.  Use to determine byte order on machine.\n"
"4 bytes integer N.  The number of filters used.\n"
"4 bytes float FLO.  The start frequency of the chirps.\n"
"4 bytes float srate.  The sample rate.\n\n"
"2xNx4 bytes float. The masses (m1,m2) m1<=m2 of the templates.\n\n"
"    | 8 bytes double time since Jan 1, 1970 UTC\n"
"    | 4 bytes int 1 if data Gaussian, 0 otherwise\n"
"    |\n"
"M X |                 |  float: Distance in Mpc at which SNR = 1\n"
"    |                 |  float: SNR\n"
"    |                 |  float: SNR with power renormaliztion\n"
"    |                 |  float: SNR with median renormalization\n"
"    | 8 x N x 4 bytes |  int: offset of waveform end (~ coalescence time)\n"
"    |                 |  float: phi (rad)=atan2(coeff 90-phase,coeff 0-phase)\n"
"    |                 |  float: r^2 splitup statistic 5.18.9 in GRASP 1.4 man\n"
"    |                 |  int: number of threshold crossings\n"
"\n";
