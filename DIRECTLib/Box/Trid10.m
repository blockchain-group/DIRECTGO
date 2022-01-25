function y = Trid10(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Trid.m
%
% Original source:
%  - http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page2904.htm
%
% Globally optimal solution:
%   f = -n*(n + 4)*(n - 1)/6
%   x(i) = [i*(10 + 1 - i)], i = 1...n
%
% Variable bounds:
%   -100 <= x(i) <= 100, i = 1...n
%   
% Problem Properties:
%   n  = 10;
%   #g = 0;
%   #h = 0;
% -------------------------------------------------------------------------
if nargin == 0
    y.nx = 10;
    y.ng = 0;
    y.nh = 0;
    y.xl = @(i) -100;
    y.xu = @(i) +100;
    y.fmin = @(i) -i*(i + 4)*(i - 1)/6;
    y.xmin = @(i) i*(10 + 1 - i);
    return
end
n = length(x);
s1 = 0;
s2 = 0;
for j = 1:n
    s1 = s1 + (x(j) - 1)^2;
end
for j = 2:n
    s2 = s2 + x(j)*x(j - 1);
end
y = s1 - s2;
end