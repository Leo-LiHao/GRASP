function [Gridout, out] = area()
% AREA
%	A Matlab version of the example: area in GRASP.
%	Example run: [Gridout, out] = area;
%
% Steve Drasco
% Summer 1998

Grid = scope_structure;

Grid.m_mn = 0.8;
Grid.m_mx = 50.0;
Grid.f_start = 140.0;

[Gridout, out] = mxTemplate_area(Grid);
