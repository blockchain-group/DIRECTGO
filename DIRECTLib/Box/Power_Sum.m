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
%   x = [1; 3; 2; 2]
%
% Variable bounds:
%   0 <= x(i) <= 4, i = 1...n
%   
% Problem Properties:
%   n  = 4;
%   #g = 0;
%   #h = 0;
% -------------------------------------------------------------------------
if nargin == 0
    y.nx = 4;
    y.ng = 0;
    y.nh = 0;
    y.xl = @(i) 0;
    y.xu = @(i) +4;
    y.fmin = @(i) 0;
    xmin = [1; 3; 2; 2];
    y.xmin = @(i) xmin(i);
    return
end
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