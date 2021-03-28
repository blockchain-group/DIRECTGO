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
%   x = zeros(n, 1);
%
% Variable bounds:
%   -15 <= x(i) <= 35, i = 1...n
%   bounds = ones(n, 1).*[-15, 35];
%
% Problem Properties:
%   n  = any dimension;
%   #g = 0;
%   #h = 0;
% -------------------------------------------------------------------------
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
y = -a*exp(-b*sqrt(1/n*s1)) - exp(1/n*s2) + a + exp(1);
end