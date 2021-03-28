function y = P6(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   P6.m
%
% Original source: 
% - Birgin, E. G., Floudas, C. A., Martínez, J. M. :Global minimization usingan 
%   Augmented Lagrangian method with variable lower-level constraints. Math. 
%   Program. Ser. A 125(1), 139–162 (2010)
%
% Globally optimal solution:
%   f* = 376.29
%   x* = (8.17001825676862, 7.56074419876632) 
%
% Constraints (including variable bounds):
%   g(1): -x(1) + ((0.2458*x(1)^2)/x(2)) + 6 <= 0;
%         0       <= x(1) <= 115.8;
%         10^(-5) <= x(2) <= 30;
%   
% Problem Properties:
%   n  = 2;
%   #g = 1;
%   #h = 0;  
% ------------------------------------------------------------------------------ 
y = 29.4*x(1)+18*x(2); 
end