function y = Csendes(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Csendes.m
%
% Original source:
%  - http://infinity77.net/global_optimization/test_functions_nd_C.html
%
% Globally optimal solution:
%   f = 0
%   x(i) = [10^(-100)], i = 1...n
%
% Variable bounds:
%   -10 <= x(i) <= 20, i = 1...n
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
    y.xu = @(i) +20;
    y.fmin = @(i) 0;
    y.xmin = @(i) 10^(-100);
    return
end
n = length(x);
y = 0;
for i = 1:n
    y = y + (x(i)^6)*(2 + sin(1/x(i)));
end
end