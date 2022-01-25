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
%   f = -2.8081311800070050^n
%   x(i) = [7.9170526915515411], i = 1...n
%
% Variable bounds:
%   0 <= x(i) <= 10, i = 1...n
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
    y.xl = @(i) 0;
    y.xu = @(i) +10;
    y.fmin = @(i) -2.8081311800070050^i;
    y.xmin = @(i) 7.9170526915515411;
    return
end
y = -prod(sqrt(x).* sin(x), 1);
end 