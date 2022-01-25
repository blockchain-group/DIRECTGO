function [y] = Perm(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Perm.m
%
% Original source:
%  - https://www.sfu.ca/~ssurjano/permdb.html
%
% Globally optimal solution:
%   f = 0
%   x = (1:n)
%
% Variable bounds:
%   -i <= x(i) <= i, i = 1...n
%   
% Problem Properties:
%   n  = any dimension;
%   #g = 0;
%   #h = 0;
% -------------------------------------------------------------------------
if nargin == 0
    y.nx = 0;
    y.ng = 0;
    y.nh = 0;
    y.xl = @(i) -i;
    y.xu = @(i) +i;
    y.fmin = @(i) 0;
    y.xmin = @(i) i;
    return
end
b = 0.5;
d = length(x);
outer = 0;

for ii = 1:d
	inner = 0;
	for jj = 1:d
        inner = inner + (jj^ii + b)*((x(jj)/jj)^ii - 1);
    end
	outer = outer + inner^2;
end
y = outer;
end