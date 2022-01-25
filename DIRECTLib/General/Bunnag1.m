function y = Bunnag1(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Bunnag1.m
%
% Original source: 
% - Bunnag D. and  Sun M. (2005, December). Genetic algorithm for 
%   constrained global optimization in continuous variables. Applied 
%   Mathematics and Computation, 171(1), 604 - 636.
%
% Globally optimal solution:
%   f* = 0.11111111111111116
%   x* = (12/9, 7/9, 4/9) 
%
% Constraints (including variable bounds):
%   g(1): x(1)-4*x(2) <= 1;
%         0 <= x(1) <= 3;
%         0 <= x(2) <= 3;
%         0 <= x(3) <= 3;
%   
% Problem Properties:
%   n  = 3;
%   #g = 1;
%   #h = 0;  
% -------------------------------------------------------------------------
if nargin == 0
    y.nx = 3;
    y.ng = 1;
    y.nh = 0;
    y.xl = @(i) 0;
    y.xu = @(i) 3;
    y.fmin = @(i) 0.11111111111111116;
    xmin = [12/9, 7/9, 4/9];
    y.xmin = @(i) xmin(i);
    y.confun = @(i) Bunnag1c(i);
    return
end
y = 9 - 8*x(1) - 6*x(2) - 4*x(3) + 2*x(1)^2 + 2*x(2)^2 + x(3)^2 +...
    2*x(1)*x(2) + 2*x(1)*x(3);
end

function [c, ceq] = Bunnag1c( x )
c   = x(1) + x(2) + 2*x(3) - 3; 
ceq = [];
end