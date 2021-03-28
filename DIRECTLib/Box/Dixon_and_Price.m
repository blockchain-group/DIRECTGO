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
%   bounds = ones(n, 1).*[-10, 10];
%   
% Problem Properties:
%   n  = any dimension;
%   #g = 0;
%   #h = 0;
% -------------------------------------------------------------------------
n = length(x);
s1 = 0;
for j = 2:n
    s1 = s1 + j*(2*x(j)^2 - x(j - 1))^2;
end
y = s1 + (x(1) - 1)^2;
end