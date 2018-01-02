function fgetoutput=oneFget()
% ONEFGET
%	This Program uses the GRASP function fget_ch() once and returns the output
%
% Steve Drasco
% Summer 1998

fgetinput.npoint=296000;
fgetinput.inlock=1;
fgetinput.seek=0;
fgetinput.calibrate=1;

%fgetinput.nchan=1;
%fgetinput.chnames={'IFO_DMRO'}

fgetinput.chnames={'IFO_DMRO','IFO_DCDM'};
fgetinput.nchan=2;

fgetoutput=mxFget_ch(fgetinput);
