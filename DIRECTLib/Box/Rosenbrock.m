function y = Rosenbrock(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Rosenbrock.m
%
% Original source:
%  - http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page2537.htm
%
% Globally optimal solution:
%   f = 0
%   x(i) = [1], i = 1...n
%   x = ones(n, 1);
%
% Variable bounds:
%   -5 <= x(i) <= 10, i = 1...n
%   bounds = ones(n, 1).*[-5, 10];
%   
% Problem Properties:
%   n  = any dimension;
%   #g = 0;
%   #h = 0;
% ------------------------------------------------------------------------------
n = length(x);
sum = 0;
for j = 1:n-1
    sum = sum + 100*(x(j)^2 - x(j + 1))^2 + (x(j) - 1)^2;
end
y = sum;
end