function y = Zakharov(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Zakharov.m
%
% Original source:
%  - http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page3088.htm
%
% Globally optimal solution:
%   f = 0
%   x(i) = [0], i = 1...n
%
% Variable bounds:
%   -5 <= x(i) <= 11, i = 1...n
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
    y.xu = @(i) +11;
    y.fmin = @(i) 0;
    y.xmin = @(i) 0;
    return
end
n = length(x);
s1 = 0;
s2 = 0;
for j = 1:n
    s1 = s1 + x(j)^2;
    s2 = s2 + 0.5*j*x(j);
end
y = s1 + s2^2 + s2^4;
end