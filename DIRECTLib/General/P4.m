function y = P4(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   P4.m
%
% Original source: 
% - Birgin, E. G., Floudas, C. A., Martínez, J. M. :Global minimization usingan 
%   Augmented Lagrangian method with variable lower-level constraints. Math. 
%   Program. Ser. A 125(1), 139–162 (2010)
%
% Globally optimal solution:
%   f* = -6.66666666
%   x* = (6, 0.66666666) 
%
% Constraints (including variable bounds):
%   g(1): x(1)*x(2)-4 <= 0;
%         0 <= x(1) <= 6;
%         0 <= x(2) <= 4;
%   
% Problem Properties:
%   n  = 2;
%   #g = 1;
%   #h = 0;  
% ------------------------------------------------------------------------------ 
y = -x(1)-x(2); 
end