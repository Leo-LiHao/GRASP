function fgetoutput=readfri()
% READFRI
%	This program reads the calibaration data from a frame file.
%	Example run: fgetoutput=readfri
%
% Steve Drasco
% Summer 1998

fgetinput.npoint=65536;
fgetinput.inlock=1;
fgetinput.seek=0;
fgetinput.calibrate=1;

fgetinput.nchan=1;
fgetinput.chnames={'IFO_DMRO'};

%fgetinput.chnames={'IFO_DMRO','IFO_DCDM'};
%fgetinput.nchan=2;

fgetoutput=mxFget_ch(fgetinput);
