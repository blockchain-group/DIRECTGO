function y = Power_Sum(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Power_Sum.m
%
% Original source:
%  - http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page670.htm
%
% Globally optimal solution:
%   f = 0
%   x = (1:n)'
%
% Variable bounds:
%   0 <= x(i) <= 4, i = 1...n
%   bounds = ones(n, 1).*[0, 4];
%   
% Problem Properties:
%   n  = any dimension;
%   #g = 0;
%   #h = 0;
% -------------------------------------------------------------------------
n = length(x);
b = [8, 18, 44, 114];
s_out = 0;
for k = 1:n
    s_in = 0;
    for j = 1:n
        s_in = s_in + x(j)^k;
    end
    s_out = s_out + (s_in - b(k))^2;
end
y = s_out;
end