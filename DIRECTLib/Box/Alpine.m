function y = Alpine(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Alpine.m
%
% Original source:
%  - http://benchmarkfcns.xyz/benchmarkfcns/alpinen2fcn.html
%
% Globally optimal solution:
%   f = -2.8081311800021023^n
%   x(i) = [7.917], i = 1...n
%   x = ones(n, 1)*7.917;
%
% Variable bounds:
%   0 <= x(i) <= 10, i = 1...n
%   bounds = ones(n, 1).*[0, 10];
%
% Problem Properties:
%   n  = any dimension;
%   #g = 0;
%   #h = 0;
% -------------------------------------------------------------------------
y = -prod(sqrt(x).* sin(x), 1);
end 