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
%   bounds = ones(2, 1).*[-4.5, 4.5];
%   
% Problem Properties:
%   n  = 2;
%   #g = 0;
%   #h = 0;
% -------------------------------------------------------------------------
y = (1.5 - x(1)*(1 - x(2)))^2 + (2.25 - x(1)*(1 - x(2)^2))^2 + (2.625 -...
    x(1)*(1 - x(2)^3))^2;
end