function fgetoutput=twoFget()
% TWOFGET
%	This program uses the GRASP function fget_ch() twice and resurns the output.
%
% Steve Drasco
% Summer 1998

fgetinput.npoint = 296000;
fgetinput.inlock = 1;
fgetinput.seek = 0;
fgetinput.calibrate = 1;

% uncomment these to get one channel
fgetinput.nchan=1;
fgetinput.chnames={'IFO_DMRO'}

% uncomment these to get two channels
%fgetinput.chnames={'IFO_DMRO','IFO_DCDM'};
%fgetinput.nchan=2;

for i = 1:2
	fgetoutput = mxFget_ch(fgetinput);
end
