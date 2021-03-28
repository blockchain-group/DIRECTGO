function y = Bunnag2(x)
% ------------------------------------------------------------------------------
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
%   f* = -6.40520658
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
% ------------------------------------------------------------------------------ 
y = x(1)^0.6+2*x(2)^0.6-2*x(2)+2*x(3)-x(4);
end