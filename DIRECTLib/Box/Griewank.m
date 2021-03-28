function y = Griewank(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Griewank.m
%
% Original source:
%  - http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page1905.htm
%
% Globally optimal solution:
%   f = 0
%   x(i) = [0], i = 1...n
%   x = zeros(n, 1);
%
% Variable bounds:
%   -600 <= x(i) <= 700, i = 1...n
%   bounds = ones(n, 1).*[-600, 700];
%   
% Problem Properties:
%   n  = any dimension;
%   #g = 0;
%   #h = 0;
% ------------------------------------------------------------------------------
n = length(x);
fr = 4000;
s = 0; 
p = 1;
for j = 1:n
    s = s + x(j)^2;
end
for j = 1:n
    p = p*cos(x(j)/sqrt(j));
end
y = s/fr - p + 1;
end