function y = Qing(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Qing.m
%
% Original source:
%  - http://infinity77.net/global_optimization/test_functions_nd_Q.html
%
% Globally optimal solution:
%   f = 0
%   x(i) = [sqrt(i)], i = 1...n
%   x = sqrt((1:10)');
%
% Variable bounds:
%   -500 <= x(i) <= 500, i = 1...n
%   bounds = ones(n, 1).*[-500, 500];
%   
% Problem Properties:
%   n  = any dimension;
%   #g = 0;
%   #h = 0;
% -------------------------------------------------------------------------
n = length(x);
y = 0;
for i = 1:n
    y = y + (x(i)^2 - i)^2;
end
end