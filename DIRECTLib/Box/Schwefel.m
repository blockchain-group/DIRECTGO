function y = Schwefel(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Schwefel.m
%
% Original source:
%  - http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page2530.htm
%
% Globally optimal solution:
%   f = 0
%   x(i) = [420.9687474737558], i = 1...n
%
% Variable bounds:
%   -500 <= x(i) <= 500, i = 1...n
%   
% Problem Properties:
%   n  = any dimension;
%   #g = 0;
%   #h = 0;
% -------------------------------------------------------------------------
if nargin == 0
    y.nx = 0;
    y.ng = 0;
    y.nh = 0;
    y.xl = @(i) -500;
    y.xu = @(i) +500;
    y.fmin = @(i) 0;
    y.xmin = @(i) 420.9687474737558;
    return
end
n = length(x);
s = sum(-x.*sin(sqrt(abs(x))));
y = 418.9828872724336*n + s;
end