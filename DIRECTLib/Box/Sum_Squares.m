function y = Sum_Squares(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Sum_Squares.m
%
% Original source:
%  - http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page674.htm
%
% Globally optimal solution:
%   f = 0
%   x(i) = [0], i = 1...n
%   x = zeros(n, 1);
%
% Variable bounds:
%   -10 <= x(i) <= 15, i = 1...n
%   bounds = ones(n, 1).*[-10, 15];
%   
% Problem Properties:
%   n  = any dimension;
%   #g = 0;
%   #h = 0;
% -------------------------------------------------------------------------
n = 2;
s = 0;
for j = 1:n
    s=s + j*x(j)^2;
end
y = s;
end