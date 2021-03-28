function y = Bukin6(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Bukin6.m
%
% Original source:
%  - https://www.sfu.ca/~ssurjano/bukin6.html
%
% Globally optimal solution:
%   f = 0
%   x = [-10; 1]
%
% Variable bounds:
%   -15 <= x(1) <= 5
%   -3  <= x(1) <= 3
%   bounds = [-15 ,5; -3, 3];
%
% Problem Properties:
%   n  = any dimension;
%   #g = 0;
%   #h = 0;
% -------------------------------------------------------------------------
term1 = 100*sqrt(abs(x(2) - 0.01*x(1)^2));
term2 = 0.01*abs(x(1) + 10);
y = term1 + term2;
end