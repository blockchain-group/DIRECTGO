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
%   f = 0.397887357729739
%   x = [pi; 2.275]
%
% Variable bounds:
%   -5 <= x(1) <= 10;
%    0 <= x(2) <= 15;
%   bounds = [-5, 10; 0, 15];
%   
% Problem Properties:
%   n  = 2;
%   #g = 0;
%   #h = 0;
% -------------------------------------------------------------------------
y = (x(2) - (5.1/(4*pi^2))*x(1)^2 + 5*x(1)/pi - 6)^2 + 10*(1 -...
    1/(8*pi))*cos(x(1)) + 10;
end