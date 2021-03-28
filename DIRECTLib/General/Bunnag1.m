function y = Bunnag1(x)
% ------------------------------------------------------------------------------
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
%   f* = 0.11111111
%   x* = (1.33333333, 0.777777777, 0.444444444) 
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
% ------------------------------------------------------------------------------ 
y = 9-8*x(1)-6*x(2)-4*x(3)+2*x(1)^2+2*x(2)^2+x(3)^2+2*x(1)*x(2)+2*x(1)*x(3);
end