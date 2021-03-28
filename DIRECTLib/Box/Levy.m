function y = Levy(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Levy.m
%
% Original source:
%  - http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page2056.htm
%
% Globally optimal solution:
%   f = 0
%   x(i) = [1], i = 1...n
%   x = ones(n, 1);
%
% Variable bounds:
%   -10 <= x(i) <= 10, i = 1...n
%   bounds = ones(n, 1).*[-10, 10];
%   
% Problem Properties:
%   n  = any dimension;
%   #g = 0;
%   #h = 0;
% ------------------------------------------------------------------------------
n = length(x);
z = zeros(1, n);
for i = 1:n
    z(i) = 1 + (x(i) - 1)/4;
end
s = sin(pi*z(1))^2;
for i = 1:n - 1
    s = s + (z(i) - 1)^2*(1 + 10*(sin(pi*z(i) + 1))^2);
end
y = s + (z(n) - 1)^2*(1 + (sin(2*pi*z(n)))^2);
end