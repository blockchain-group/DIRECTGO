function y = Matyas(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Matyas.m
%
% Original source:
%  - http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page2213.htm
%
% Globally optimal solution:
%   f = 0
%   x = [0, 0]
%
% Variable bounds:
%   -10 <= x(i) <= 15, i = 1...2
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
    y.xl = @(i) -10;
    y.xu = @(i) +15;
    y.fmin = @(i) 0;
    y.xmin = @(i) 0;
    return
end
y = 0.26*(x(1)^2 + x(2)^2) - 0.48*x(1)*x(2);
end