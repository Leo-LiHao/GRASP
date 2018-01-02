function out=cubic_grid_structure()
% CUBIC_GRID_STRUCTURE
%	Output is an empty cubic_grid structure.
%
% Steve Drasco
% Summer 1998

out = struct('n',[],'m_mn',[],'m_mx',[],'dm',[],'match',[],'angle',[], ...
	'order',[],'srate',[],'flo',[],'ftau',[],'detector',[],'coef',[]);
