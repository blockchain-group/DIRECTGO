function y = Schwefel(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Schwefel.m
%
% Original source:
%  - http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page2530.htm
%
% Globally optimal solution:
%   f = 0
%   x(i) = [0], i = 1...n
%   x = zeros(n, 1);
%
% Variable bounds:
%   -500 <= x(i) <= 500, i = 1...n
%   bounds = ones(n, 1).*[-500, 500];
%   
% Problem Properties:
%   n  = any dimension;
%   #g = 0;
%   #h = 0;
% ------------------------------------------------------------------------------
n = length(x);
s = sum(-x.*sin(sqrt(abs(x))));
y = 418.98288727245*n + s;
end