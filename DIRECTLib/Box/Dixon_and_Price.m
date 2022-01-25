function y = Dixon_and_Price(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Dixon_and_Price.m
%
% Original source:
%  - http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page1240.htm
%
% Globally optimal solution:
%   f = 0
%   x(i) = [2^(-((2^i - 2)/(2^i)))] i = 1...n;
%
% Variable bounds:
%   -10 <= x(i) <= 10, i = 1...n
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
    y.xu = @(i) +10;
    y.fmin = @(i) 0;
    y.xmin = @(i) 2^(-((2^i - 2)/(2^i)));
    return
end
n = length(x);
s1 = 0;
for j = 2:n
    s1 = s1 + j*(2*x(j)^2 - x(j - 1))^2;
end
y = s1 + (x(1) - 1)^2;
end