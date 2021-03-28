function y = P8(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   P8.m
%
% Original source: 
% - Birgin, E. G., Floudas, C. A., Martínez, J. M. :Global minimization usingan 
%   Augmented Lagrangian method with variable lower-level constraints. Math. 
%   Program. Ser. A 125(1), 139–162 (2010)
%
% Globally optimal solution:
%   f* = -118.7
%   x* = (-3.17361715437436, 1.72471153374410) 
%
% Constraints (including variable bounds):
%   g(1): x(2)-x(1)^2-2*x(1)+2 <= 0;
%   g(2): -x(1)+x(2)-8         <= 0;
%         -8 <= x(1) <= 10;
%         -8 <= x(2) <= 10;
%   
% Problem Properties:
%   n  = 2;
%   #g = 2;
%   #h = 0;  
% ------------------------------------------------------------------------------ 
y = x(1)^4-14*x(1)^2+24*x(1)-x(2)^2;
end