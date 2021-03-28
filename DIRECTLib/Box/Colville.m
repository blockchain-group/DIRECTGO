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
%   x = ones(4, 1);
%
% Variable bounds:
%   -10 <= x(i) <= 10, i = 1...4
%   bounds = ones(4, 1).*[-10, 10];
%   
% Problem Properties:
%   n  = 4;
%   #g = 0;
%   #h = 0;
% -------------------------------------------------------------------------
y = 100*(x(1)^2 - x(2))^2 + (x(1) - 1)^2 + (x(3) - 1)^2 + 90*(x(3)^2 -...
    x(4))^2 + 10.1*((x(2) - 1)^2 + (x(4) - 1)^2) + 19.8*(x(2) -...
    1)*(x(4) - 1);
end    