function y = Levy(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Levy.m
%
% Original source:
%  - http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page2056.htm
%
% Globally optimal solution:
%   f = 0
%   x = [1], i = 1...n
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
    y.xmin = @(i) 1;
    return
end
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