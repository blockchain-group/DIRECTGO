function y = Beale(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Beale.m
%
% Original source:
%  - http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page288.htm
%
% Globally optimal solution:
%   f = 0
%   x = [3; 0.5]
%
% Variable bounds:
%   -4.5 <= x(i) <= 4.5, i = 1...n
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
    y.xl = @(i) -4.5;
    y.xu = @(i) +4.5;
    y.fmin = @(i) 0;
    xmin = [3; 0.5];
    y.xmin = @(i) xmin(i);
    return
end
y = (1.5 - x(1)*(1 - x(2)))^2 + (2.25 - x(1)*(1 - x(2)^2))^2 + (2.625 -...
    x(1)*(1 - x(2)^3))^2;
end