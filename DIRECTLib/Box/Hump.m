function y = Hump(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Hump.m
%
% Original source:
%  - http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page1621.htm
%
% Globally optimal solution:
%   f = 0
%   x = [0.0898; -0.7126]
%
% Variable bounds:
%   -5 <= x(i) <= 5, i = 1...2
%   bounds = ones(2, 1).*[-5, 5];
%   
% Problem Properties:
%   n  = 2;
%   #g = 0;
%   #h = 0;
% -------------------------------------------------------------------------
y = 4*x(1)^2 - 2.1*x(1)^4 + x(1)^6/3 + x(1)*x(2) - 4*x(2)^2 + 4*x(2)^4;
end