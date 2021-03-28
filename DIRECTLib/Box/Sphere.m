function y = Sphere(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Sphere.m
%
% Original source:
%  - http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page1113.htm
%
% Globally optimal solution:
%   f = 0
%   x(i) = [0], i = 1...n
%   x = zeros(n, 1);
%
% Variable bounds:
%   -5.12 <= x(i) <= 6.12, i = 1...n
%   bounds = ones(n, 1).*[-5.12, 6.12];
%   
% Problem Properties:
%   n  = any dimension;
%   #g = 0;
%   #h = 0;
% -------------------------------------------------------------------------
n = length(x);
s = 0;
for j = 1:n
    s = s + x(j)^2;
end
y = s;
end