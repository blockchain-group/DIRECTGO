function y = Rotated_Hyper_Ellipsoid(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Rotated_Hyper_Ellipsoid.m
%
% Original source:
%  - https://www.sfu.ca/~ssurjano/rothyp.html
%
% Globally optimal solution:
%   f = 0
%   x(i) = [0], i = 1...n
%   x = zeros(n, 1);
%
% Variable bounds:
%   -65.536 <= x(i) <= 66.536, i = 1...n
%   bounds = ones(n, 1).*[-65.536, 66.536];
%    
% Problem Properties:
%   n  = any dimension;
%   #g = 0;
%   #h = 0;
% -------------------------------------------------------------------------
d = length(x);
outer = 0;
for ii = 1:d
    inner = 0;
    for jj = 1:ii
        xj = x(jj);
        inner = inner + xj^2;
    end
    outer = outer + inner;
end
y = outer;
end