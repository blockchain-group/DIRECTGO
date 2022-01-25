function y = Bohachecsky1(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Bohachecsky1.m
%
% Original source:
%  - http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page595.htm
%
% Globally optimal solution:
%   f = 0
%   x = [0; 0]
%
% Variable bounds:
%   -100 <= x(i) <= 110, i = 1...2
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
    y.xl = @(i) -100;
    y.xu = @(i) +110;
    y.fmin = @(i) 0;
    y.xmin = @(i) 0;
    return
end
y = x(1)^2 + 2*x(2)^2 - 0.3*cos(3*pi*x(1)) - 0.4*cos(4*pi*x(2)) + 0.7;
end