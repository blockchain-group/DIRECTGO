function y = Bunnag2(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Bunnag2.m
%
% Original source: 
% - Bunnag D. and  Sun M. (2005, December). Genetic algorithm for 
%   constrained global optimization in continuous variables. Applied 
%   Mathematics and Computation, 171(1), 604 - 636.
%
% Globally optimal solution:
%   f* = -6.4052065800118604954604961676523
%   x* = (1, 4, 0, 4) 
%
% Constraints (including variable bounds):
%   g(1): x(1)+2*x(3)  <= 4;
%   g(2): -3*x(1)+x(4) <= 1;
%         0 <= x(1) <= 4;
%         0 <= x(2) <= 4;
%         0 <= x(3) <= 4;
%         0 <= x(4) <= 4;
%   
% Problem Properties:
%   n  = 4;
%   #g = 2;
%   #h = 0;  
% -------------------------------------------------------------------------
if nargin == 0
    y.nx = 4;
    y.ng = 2;
    y.nh = 0;
    y.xl = @(i) 0;
    y.xu = @(i) 4;
    y.fmin = @(i) -6.4052065800118604954604961676523;
    xmin = [1, 4, 0, 4];
    y.xmin = @(i) xmin(i);
    y.confun = @(i) Bunnag2c(i);
    return
end
y = x(1)^0.6 + 2*x(2)^0.6 - 2*x(2) + 2*x(3) - x(4);
end

function [c, ceq] = Bunnag2c( x )
c(1) = x(1) + 2*x(3) - 4; 
c(2) = -3*x(1) + x(4) - 1; 
ceq = [];
end