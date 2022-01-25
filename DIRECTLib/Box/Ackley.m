function y = Ackley(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Ackley.m
%
% Original source:
%  - http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page295.htm
%
% Globally optimal solution:
%   f = 0
%   x(i) = [0], i = 1...n
%
% Variable bounds:
%   -15 <= x(i) <= 35, i = 1...n
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
    y.xl = @(i) -15;
    y.xu = @(i) +35;
    y.fmin = @(i) 0;
    y.xmin = @(i) 0;
    return
end
n = length(x);
a = 20; 
b = 0.2; 
c = 2*pi;
s1 = 0; 
s2 = 0;
for i=1:n
	s1 = s1 + x(i)^2;
	s2 = s2 + cos(c*x(i));
end
y = -a*exp(-b*(1/n*s1)^(1/2)) - exp(1/n*s2) + a + exp(1);
end