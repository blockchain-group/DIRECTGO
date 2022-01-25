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
    xl = [-15; -3];
    y.xl = @(i) xl(i);
    xu = [5; 3];
    y.xu = @(i) xu(i);
    y.fmin = @(i) 0;
    xmin = [-10; 1];
    y.xmin = @(i) xmin(i);
    return
end
term1 = 100*sqrt(abs(x(2) - 0.01*x(1)^2));
term2 = 0.01*abs(x(1) + 10);
y = term1 + term2;
end