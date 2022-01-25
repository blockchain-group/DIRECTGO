function y = Langermann(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Langermann.m
%
% Original source:
%  - https://www.sfu.ca/~ssurjano/langer.html
%
% Globally optimal solution:
%   f = -4.1558092918
%   x = [2.79340196434474; 1.59723280665210]
%
% Variable bounds:
%   0 <= x(i) <= 10, i = 1...2
%
% Problem Properties:
%   n  = 2;
%   #g = 0;
%   #h = 0;
% -------------------------------------------------------------------------
if nargin == 0
    y.nx = 2;
    y.ng = 0;
    y.nh = 0;
    y.xl = @(i) 0;
    y.xu = @(i) +10;
    y.fmin = @(i) -4.155809291843469;
    xmin = [2.79340196434474;1.59723280665210];
    y.xmin = @(i) xmin(i);
    return
end
d = length(x);
m = 5;
c = [1, 2, 5, 2, 3];
A = [3, 5; 5, 2; 2, 1; 1, 4; 7, 9];
outer = 0;
for ii = 1:m
    inner = 0;
    for jj = 1:d
        xj = x(jj);
        Aij = A(ii, jj);
        inner = inner + (xj - Aij)^2;
    end
    new = c(ii) * exp(-inner/pi) * cos(pi*inner);
    outer = outer + new;
end
y = outer;
end
