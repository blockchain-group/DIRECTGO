function y = Sum_Square(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Sum_Square.m
%
% Original source:
%  - https://www.sfu.ca/~ssurjano/sumsqu.html
%
% Globally optimal solution:
%   f = 0
%   x(i) = [0], i = 1...n
%
% Variable bounds:
%   -10 <= x(i) <= 15, i = 1...n
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
    y.xl = @(i) -10;
    y.xu = @(i) +15;
    y.fmin = @(i) 0;
    y.xmin = @(i) 0;
    return
end
n = length(x);
s = 0;
for j = 1:n
    s = s + j*x(j)^2;
end
y = s;
return