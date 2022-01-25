function y = Booth(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Booth.m
%
% Original source:
%  - http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page816.htm
%
% Globally optimal solution:
%   f = 0
%   x = [1; 3]
%
% Variable bounds:
%   -10 <= x(i) <= 10, i = 1...2
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
    y.xu = @(i) +10;
    y.fmin = @(i) 0;
    xmin = [1; 3];
    y.xmin = @(i) xmin(i);
    return
end
y = (x(1) + 2*x(2) - 7)^2 + (2*x(1) + x(2) - 5)^2;
end