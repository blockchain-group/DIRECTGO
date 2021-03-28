function y = P16(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   P16.m
%
% Original source: 
% - Birgin, E. G., Floudas, C. A., Martínez, J. M. :Global minimization usingan 
%   Augmented Lagrangian method with variable lower-level constraints. Math. 
%   Program. Ser. A 125(1), 139–162 (2010)
%
% Globally optimal solution:
%   f* = 0.704902358236
%   x* = (1.820175974329, 2.956011468603) 
%
% Constraints (including variable bounds):
%   g(1): ((x(1)-1)/(36-12*x(1)))-1.5834  <= 0;
%   g(2): ((x(2)-x(1))/(32-8*x(2)))-3.625 <= 0;
%   g(3): ((5-x(2))/4)-1                  <= 0;
%   g(4): -((x(1)-1)/(36-12*x(1)))        <= 0;
%   g(5): -((x(2)-x(1))/(32-8*x(2)))      <= 0;
%   g(6): -((5-x(2))/4)                   <= 0;
%         1 <= x(1) <= 3;
%         1 <= x(2) <= 4;
%   
% Problem Properties:
%   n  = 2;
%   #g = 6;
%   #h = 0;  
% ------------------------------------------------------------------------------ 
y = ((x(1)-1)/(36-12*x(1)))+((x(2)-x(1))/(32-8*x(2)))+((5-x(2))/4); 
end