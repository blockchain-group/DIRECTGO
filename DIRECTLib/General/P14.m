function y = P14(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   P14.m
%
% Original source: 
% - Birgin, E. G., Floudas, C. A., Martínez, J. M. :Global minimization usingan 
%   Augmented Lagrangian method with variable lower-level constraints. Math. 
%   Program. Ser. A 125(1), 139–162 (2010)
%
% Globally optimal solution:
%   f* = -4.5143835194
%   x* = (1.333333333333, 4, 0) 
%
% Constraints (including variable bounds):
%   g(1): (1/3)*x(2)-x(1)-2          <= 0;
%   g(2): x(1)+2*((1/3)*x(2)-x(1))-4 <= 0;
%   g(3): x(2)+2*x(3)-4              <= 0;
%   g(4): -((1/3)*x(2)-x(1))         <= 0;
%         10^(-5) <= x(1) <= 3;
%         10^(-5) <= x(2) <= 4;
%         0       <= x(3) <= 1;
%   
% Problem Properties:
%   n  = 3;
%   #g = 4;
%   #h = 0;  
% ------------------------------------------------------------------------------ 
y = x(1)^(0.6)+x(2)^(0.6)-2*x(1)-(4/3)*x(2)+3*x(3); 
end