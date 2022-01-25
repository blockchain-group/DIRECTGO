function y = Sum_of_Different_Powers(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Sum_of_Different_Powers.m
%
% Original source:
%  - https://www.sfu.ca/~ssurjano/sumpow.html
%
% Globally optimal solution:
%   f = 0
%   x(i) = [0], i = 1...n
%
% Variable bounds:
%   -1 <= x(i) <= 2.5, i = 1...n
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
    y.xl = @(i) -1;
    y.xu = @(i) +2.5;
    y.fmin = @(i) 0;
    y.xmin = @(i) 0;
    return
end
d = length(x);
sum = 0;
for ii = 1:d
    new = (abs(x(ii)))^(ii + 1);
    sum = sum + new;
end
y = sum;
end