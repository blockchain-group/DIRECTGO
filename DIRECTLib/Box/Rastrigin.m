function y = Rastrigin(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Rastrigin.m
%
% Original source:
%  - http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page2607.htm
%
% Globally optimal solution:
%   f = 0
%   x(i) = [0], i = 1...n
%   x = zeros(n, 1);
%
% Variable bounds:
%   -6.12 <= x(i) <= 5.12, i = 1...n
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
    y.xl = @(i) -6.12;
    y.xu = @(i) +5.12;
    y.fmin = @(i) 0;
    y.xmin = @(i) 0;
    return
end
n = length(x);
s = 0;
for j = 1:n
    s = s + (x(j)^2 - 10*cos(2*pi*x(j)));
end
y = 10*n + s;
end