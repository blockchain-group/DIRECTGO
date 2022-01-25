function y = Styblinski_Tang(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Styblinski_Tang.m
%
% Original source:
%  - https://www.sfu.ca/~ssurjano/stybtang.html
%
% Globally optimal solution:
%   f = -39.166165703771426*n
%   x(i) = [-2.9035340311065125], i = 1...n
%
% Variable bounds:
%   -5 <= x(i) <= 5, i = 1...n
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
    y.xl = @(i) -5;
    y.xu = @(i) +5;
    y.fmin = @(i) -39.166165703771426*i;
    y.xmin = @(i) -2.9035340311065125;
    return
end
d = length(x);
sum = 0;
for ii = 1:d
	xi = x(ii);
	new = xi^4 - 16*xi^2 + 5*xi;
	sum = sum + new;
end
y = sum/2;
end