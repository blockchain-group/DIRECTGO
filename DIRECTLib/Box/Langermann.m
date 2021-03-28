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
%   x(i) = [0], i = 1...n
%   x = zeros(n, 1);
%
% Variable bounds:
%   0 <= x(i) <= 10, i = 1...2
%   bounds = ones(2, 1).*[0, 10];
%
% Problem Properties:
%   n  = 2;
%   #g = 0;
%   #h = 0;
% -------------------------------------------------------------------------
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
