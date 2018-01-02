function [fgetoutput, response]=power_spectrumF()
% POWER_SPECTRUMF
%	A Matlab version of the GRASP example: power_spectrumF.
%	Example run: [fgetoutput, response]=power_spectrumF;
%
% Steve Drasco
% Summer 1998

npoint = 65536;
fgetinput.nchan = 1;
fgetinput.chnames = {'IFO_DMRO'};
fgetinput.inlock = 1;
fgetinput.npoint = npoint;
fgetinput.calibrate = 1;

% I won't skip any data
fgetinput.seek = 0;
fgetoutput = mxFget_ch(fgetinput);
srate = fgetoutput.srate;
data = fft(fgetoutput.data{1});
response = mxGRnormalize(fgetoutput.fri,fgetoutput.frinum, npoint, srate);

% one-sided power-spectrum normalization, to get meters/rHz
factor = sqrt(2/(srate*npoint));

% frequency
freq = (1:(npoint/2)-1)*srate/npoint;

% real and imaginary parts of tilde c0
c0_real = real( data( 2:npoint/2 ) );
c0_imag = imag( data( 2:npoint/2 ) );

% real and imaginary parts of R
res_real = response(   1 + 2*( 2:npoint/2  )   );
res_imag = response(   2 + 2*( 2:npoint/2  )   );

% real and imaginary parts of tilde dl
dl_real = c0_real.*res_real - c0_imag.*res_imag;
dl_imag = c0_real.*res_imag + c0_imag.*res_real;
spectrum = factor * sqrt(dl_real.^2 + dl_imag.^2);

% plot the results
loglog(freq, spectrum);
xlabel('frequency (Hz)');
ylabel('noise m/(Hz^1/2)'); 
