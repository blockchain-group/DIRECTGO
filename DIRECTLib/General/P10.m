function y = P10(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   P10.m
%
% Original source: 
% - Birgin, E. G., Floudas, C. A., Martínez, J. M. :Global minimization usingan 
%   Augmented Lagrangian method with variable lower-level constraints. Math. 
%   Program. Ser. A 125(1), 139–162 (2010)
%
% Globally optimal solution:
%   f* = 0.74178
%   x* = (0.129409522551260, 0.482962913144536) 
%
% Constraints (including variable bounds):
%   g(1): -16*x(1)*x(2)+1      <= 0;
%   g(2): -4*x(1)^2-4*x(2)^2+1 <= 0;
%         0 <= x(1) <= 1;
%         0 <= x(2) <= 1;
%   
% Problem Properties:
%   n  = 2;
%   #g = 2;
%   #h = 0;  
% ------------------------------------------------------------------------------ 
y = 2*x(1)+x(2); 
end