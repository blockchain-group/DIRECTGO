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
%   f = -39.1661657038*n
%   x(i) = [-2.903534], i = 1...n
%   x = ones(n, 1)*(-2.903534);
%
% Variable bounds:
%   -5 <= x(i) <= 5, i = 1...n
%   bounds = ones(n, 1).*[-5, 5];
%   
% Problem Properties:
%   n  = any dimension;
%   #g = 0;
%   #h = 0;
% -------------------------------------------------------------------------
d = length(x);
sum = 0;
for ii = 1:d
	xi = x(ii);
	new = xi^4 - 16*xi^2 + 5*xi;
	sum = sum + new;
end
y = sum/2;
end