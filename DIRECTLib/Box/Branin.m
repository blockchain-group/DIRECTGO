function y = Branin(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Branin.m
%
% Original source:
%  - http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page913.htm
%
% Globally optimal solution:
%   f = 0.39788735772973815
%   x = [3.14159264890551; 2.275000033046208]
%
% Variable bounds:
%   -5 <= x(1) <= 10;
%    0 <= x(2) <= 15;
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
    xl = [-5; 0];
    y.xl = @(i) xl(i);
    xu = [10; 15];
    y.xu = @(i) xu(i);
    y.fmin = @(i) 0.39788735772973815;
    xmin = [3.14159264890551; 2.275000033046208];
    y.xmin = @(i) xmin(i);
    return
end
y = (x(2) - (5.1/(4*pi^2))*x(1)^2 + 5*x(1)/pi - 6)^2 + 10*(1 -...
    1/(8*pi))*cos(x(1)) + 10;
end