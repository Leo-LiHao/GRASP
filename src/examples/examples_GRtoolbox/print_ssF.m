function print_ssF();
% PRINT_SSF
%	A Matlab version of the GRASP example: print_ssF.
%	Example run: print_ssF
%
% Steve Drasco
% Summer 1998

% get some data
fgetinput.npoint = 256;
fgetinput.nchan = 1;
fgetinput.chnames = {'IFO_DMRO'};
fgetinput.inlock = 0;
fgetinput.seek = 1;
fgetinput.calibrate = 1;
fgetoutput = mxFget_ch(fgetinput);

% call mxGRcalibrate
srate = fgetoutput.srate;
npoint = 4096;
cplx = mxGRcalibrate(fgetoutput.fri, fgetoutput.frinum, npoint, srate, 2, 0);

% plot output
freq=1:npoint/2;
freq = freq*srate/npoint;
imaginary = cplx(2*(1:npoint/2));
real = cplx(2*(1:npoint/2)+1);
plot(freq, real, 'b', freq, imaginary, 'r');
