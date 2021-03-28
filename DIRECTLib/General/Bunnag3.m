function y = Bunnag3(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Bunnag3.m
%
% Original source: 
% - Bunnag D. and  Sun M. (2005, December). Genetic algorithm for 
%   constrained global optimization in continuous variables. Applied 
%   Mathematics and Computation, 171(1), 604 - 636.
%
% Globally optimal solution:
%   f* = -16.36904353
%   x* = (0, 0, 4, 1.33322012, 0)
%
% Constraints (including variable bounds):
%   g(1): x(1)+2*x(4)          <= 4;
%   g(2): 3*x(1)+3*x(4)+x(5)   <= 4;
%   g(3): 2*x(2)+4*x(4)+2*x(5) <= 6;
%         0 <= x(1) <= 3;
%         0 <= x(2) <= 2;
%         0 <= x(3) <= 4;
%         0 <= x(4) <= 4;
%         0 <= x(5) <= 2;
%   
% Problem Properties:
%   n  = 5;
%   #g = 3;
%   #h = 0;  
% ------------------------------------------------------------------------------ 
y = x(1)^0.6+x(2)^0.6+x(3)^0.6-4*x(3)-2*x(4)+5*x(5);
end