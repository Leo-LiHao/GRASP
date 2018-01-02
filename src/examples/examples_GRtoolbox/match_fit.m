function [semimajor, semiminor, theta, mcoef, tstp]=match_fit(m1, m2, matchcont, order)
% MATCH_FIT
%	A Matlab version of the GRASP example: match_fit.
%	Example run: [semimajor, semiminor, theta, mcoef, tstp]=match_fit(m1, m2, matchcont, order)
%
% Steve Drasco
% Summer 1998

srate = 50000;
detector_num = 15;
flo = 120;
ftau = 140;

[site_parameters, site_name, noise_file, whiten_file]=mxDetector_site('detectors.dat', detector_num);

[semimajor, semiminor, theta, mcoef, tstp] = mxMatch_parab(m1, m2, matchcont, order, srate, flo, ftau, noise_file);
if tstp
	semimajor
	semiminor
	theta
	mcoef
elseif ~tstp
	[semimajor, semiminor, theta, mcoef, tstp] = mxMatch_parab(m1, m2, matchcont, order, srate, flo, ftau, noise_file);
	semimajor
        semiminor
        theta
        mcoef
end
