function y = Holder(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Holder.m
%
% Original source:
%  - https://www.sfu.ca/~ssurjano/holder.html
%
% Globally optimal solution:
%   f = -19.2085
%   x(i) = [8.05502, 9.66459]
%
% Variable bounds:
%   -10 <= x(i) <= 10, i = 1...2
%   bounds = ones(2, 1).*[-10, 10];
%
% Problem Properties:
%   n  = 2;
%   #g = 0;
%   #h = 0;
% -------------------------------------------------------------------------
fact1 = sin(x(1))*cos(x(2));
fact2 = exp(abs(1 - sqrt(x(1)^2 + x(2)^2)/pi));
y = -abs(fact1*fact2);
end