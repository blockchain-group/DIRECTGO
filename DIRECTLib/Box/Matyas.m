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
%   bounds = ones(2, 1).*[-10, 15];
%   
% Problem Properties:
%   n  = 2;
%   #g = 0;
%   #h = 0;
% -------------------------------------------------------------------------
y = 0.26*(x(1)^2 + x(2)^2) - 0.48*x(1)*x(2);
end