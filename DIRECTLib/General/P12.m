function y = P12(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   P12.m
%
% Original source: 
% - Birgin, E. G., Floudas, C. A., Martínez, J. M. :Global minimization usingan 
%   Augmented Lagrangian method with variable lower-level constraints. Math. 
%   Program. Ser. A 125(1), 139–162 (2010)
%
% Globally optimal solution:
%   f* = -16.7392985937000
%   x* = (0.717536188588019) 
%
% Constraints (including variable bounds):
%   g(1): 2-2*x(1)^4-3  <= 0;
%   g(2): -(2-2*x(1)^4) <= 0;
%         0 <= x(1) <= 2;
%   
% Problem Properties:
%   n  = 1;
%   #g = 2;
%   #h = 0;  
% ------------------------------------------------------------------------------ 
y = -12*x(1)+6*x(1)^4+4*x(1)^8-10; 
end