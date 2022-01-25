function y = Bunnag7(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Bunnag7.m
%
% Original source: 
% - Bunnag D. and  Sun M. (2005, December). Genetic algorithm for 
%   constrained global optimization in continuous variables. Applied 
%   Mathematics and Computation, 171(1), 604 - 636.
%
% Globally optimal solution:
%   f* = -39
%   x* = (1, 0 , 0, 1, 1, 1, 0, 1, 1, 1)
%
% Constraints (including variable bounds):
%   g(1): -2*x(1) - 6*x(2) - x(3) - 3*x(5) - 3*x(6) - 2*x(7) - 6*x(8) - 2*x(9) - 2*x(10) <= -4;
%   g(2): 6*x(1) - 5*x(2) + 8*x(3) - 3*x(4) + x(6) + 3*x(7) + 8*x(8) + 9*x(9) - 3*x(10) <= 22;
%   g(3): -5*x(1) + 6*x(2) + 5*x(3) + 3*x(4) + 8*x(5) - 8*x(6) + 9*x(7) + 2*x(8) - 9*x(10) <= -6; 
%   g(4): 9*x(1) + 5*x(2) - 9*x(4) + x(5) - 8*x(6) + 3*x(7) - 9*x(8) - 9*x(9) - 3*x(10) <= -23;
%   g(5): -8*x(1) + 7*x(2) - 4*x(3) - 5*x(4) - 9*x(5) + x(6) - 7*x(7) - x(8) + 3*x(9) - 2*x(10) <= -12;
%         0 <= x(1) <= 1;
%         0 <= x(2) <= 1;
%         0 <= x(3) <= 1;
%         0 <= x(4) <= 1;
%         0 <= x(5) <= 1;
%         0 <= x(6) <= 1;
%         0 <= x(7) <= 1;
%         0 <= x(8) <= 1;
%         0 <= x(9) <= 1;
%         0 <= x(10) <= 1;
%   
% Problem Properties:
%   n  = 10;
%   #g = 5;
%   #h = 0;  
% -------------------------------------------------------------------------
if nargin == 0
    y.nx = 10;
    y.ng = 5;
    y.nh = 0;
    y.xl = @(i) 0;
    y.xu = @(i) 1;
    y.fmin = @(i) -39;
    xmin = [1, 0 , 0, 1, 1, 1, 0, 1, 1, 1];
    y.xmin = @(i) xmin(i);
    y.confun = @(i) Bunnag7c(i);
    return
end
d = [48, 42, 48, 45, 44, 41, 47, 42, 45, 46];
b = transpose(d);
y = sum(b.*x) - 50*(sum((x.^2)));  
end

function [c, ceq] = Bunnag7c( x )
c(1) = -2*x(1) - 6*x(2) - x(3) - 3*x(5) - 3*x(6) - 2*x(7) - 6*x(8) - 2*x(9) - 2*x(10) + 4;
c(2) = 6*x(1) - 5*x(2) + 8*x(3) - 3*x(4) + x(6) + 3*x(7) + 8*x(8) + 9*x(9) - 3*x(10) - 22;
c(3) = -5*x(1) + 6*x(2) + 5*x(3) + 3*x(4) + 8*x(5) - 8*x(6) + 9*x(7) + 2*x(8) - 9*x(10) + 6;
c(4) = 9*x(1) + 5*x(2) - 9*x(4) + x(5) - 8*x(6) + 3*x(7) - 9*x(8) - 9*x(9) - 3*x(10) + 23;
c(5) = -8*x(1) + 7*x(2) - 4*x(3) - 5*x(4) - 9*x(5) + x(6) - 7*x(7) - x(8) + 3*x(9) - 2*x(10) + 12;
ceq = [];
end