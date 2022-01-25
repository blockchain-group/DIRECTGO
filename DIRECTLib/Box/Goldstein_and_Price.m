function y = Goldstein_and_Price(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Goldstein_and_Price.m
%
% Original source:
%  - http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page1760.htm
%
% Globally optimal solution:
%   f = 3
%   x = [0; -1]
%
% Variable bounds:
%   -2 <= x(i) <= 2, i = 1...2
%   
% Problem Properties:
%   n  = 2;
%   #g = 0;
%   #h = 0;
% -------------------------------------------------------------------------
if nargin == 0
    y.nx = 2;
    y.ng = 0;
    y.nh = 0;
    y.xl = @(i) -2;
    y.xu = @(i) +2;
    y.fmin = @(i) 3;
    xmin = [0; -1];
    y.xmin = @(i) xmin(i);
    return
end
a = 1 + (x(1) + x(2) + 1)^2*(19 - 14*x(1) + 3*x(1)^2 - 14*x(2) +...
    6*x(1)*x(2) + 3*x(2)^2);
b = 30 + (2*x(1) - 3*x(2))^2*(18 - 32*x(1) + 12*x(1)^2 + 48*x(2) -...
    36*x(1)*x(2) + 27*x(2)^2);
y = a*b;
end