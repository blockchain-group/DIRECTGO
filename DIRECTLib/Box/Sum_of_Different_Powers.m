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
%   x = zeros(n, 1);
%
% Variable bounds:
%   -1 <= x(i) <= 2.5, i = 1...n
%   bounds = ones(n, 1).*[-1, 2.5];
%   
% Problem Properties:
%   n  = any dimension;
%   #g = 0;
%   #h = 0;
% -------------------------------------------------------------------------
d = length(x);
sum = 0;
for ii = 1:d
    new = (abs(x(ii)))^(ii + 1);
    sum = sum + new;
end
y = sum;
end