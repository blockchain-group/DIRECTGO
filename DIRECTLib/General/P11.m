function y = P11(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   P11.m
%
% Original source: 
% - Birgin, E. G., Floudas, C. A., Martínez, J. M. :Global minimization usingan 
%   Augmented Lagrangian method with variable lower-level constraints. Math. 
%   Program. Ser. A 125(1), 139–162 (2010)
%
% Globally optimal solution:
%   f* = -0.5
%   x* = (0.5, 0.5) 
%
% Constraints (including variable bounds):
%   g(1): 4*x(1)*x(2)+2*x(1)+2*x(2)-3 <= 0;
%         0 <= x(1) <= 1;
%         0 <= x(2) <= 1;
%   
% Problem Properties:
%   n  = 2;
%   #g = 1;
%   #h = 0;  
% ------------------------------------------------------------------------------ 
y = -2*x(1)*x(2); 
end