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
%   x(i) = [0], i = 1...n
%   x = zeros(n, 1);
%
% Variable bounds:
%   -10 <= x(i) <= 20, i = 1...n
%   bounds = ones(n, 1).*[-10, 20];
%
% Problem Properties:
%   n  = any dimension;
%   #g = 0;
%   #h = 0;
% -------------------------------------------------------------------------
n = length(x);
y = 0;
for i = 1:n
    y = y + (x(i)^6)*(2 + sin(1/x(i)));
end
end