function y = Permdb(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Permdb.m
%
% Original source:
%  - https://www.sfu.ca/~ssurjano/perm0db.html
%
% Globally optimal solution:
%   f = 0
%   x = ones(n, 1)./(1:n)
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
    y.xmin = @(i) 1/i;
    return
end
b = 10;
d = length(x);
outer = 0;
for ii = 1:d
	inner = 0;
	for jj = 1:d
        inner = inner + (jj + b)*(x(jj)^ii - (1/jj)^ii);
	end
	outer = outer + inner^2;
end
y = outer;
end