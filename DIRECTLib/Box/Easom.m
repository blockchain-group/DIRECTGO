function y = Easom(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Easom.m
%
% Original source:
%  - http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page1361.htm
%
% Globally optimal solution:
%   f = -1
%   x = [pi; pi]
%
% Variable bounds:
%   -100 <= x(i) <= 100, i = 1...2
%   bounds = ones(2, 1).*[-100, 100]
%   
% Problem Properties:
%   n  = 2;
%   #g = 0;
%   #h = 0;
% -------------------------------------------------------------------------
y = -cos(x(1))*cos(x(2))*exp(-(x(1) - pi)^2 - (x(2) - pi)^2);
end