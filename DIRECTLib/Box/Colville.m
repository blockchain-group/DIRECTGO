function y = Colville(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Colville.m
%
% Original source:
%  - http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page1016.htm
%
% Globally optimal solution:
%   f = 0
%   x(i) = [1], i = 1...4
%
% Variable bounds:
%   -10 <= x(i) <= 10, i = 1...4
%   
% Problem Properties:
%   n  = 4;
%   #g = 0;
%   #h = 0;
% -------------------------------------------------------------------------
if nargin == 0
    y.nx = 4;
    y.ng = 0;
    y.nh = 0;
    y.xl = @(i) -10;
    y.xu = @(i) +10;
    y.fmin = @(i) 0;
    y.xmin = @(i) 1;
    return
end
y = 100*(x(1)^2 - x(2))^2 + (x(1) - 1)^2 + (x(3) - 1)^2 + 90*(x(3)^2 -...
    x(4))^2 + 10.1*((x(2) - 1)^2 + (x(4) - 1)^2) + 19.8*(x(2) -...
    1)*(x(4) - 1);
end    