function y = P5(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   P5.m
%
% Original source: 
% - Birgin, E. G., Floudas, C. A., Martínez, J. M. :Global minimization usingan 
%   Augmented Lagrangian method with variable lower-level constraints. Math. 
%   Program. Ser. A 125(1), 139–162 (2010)
%
% Globally optimal solution:
%   f* = 201.16
%   x* = (6.29343204289132, 3.82183760411190) 
%
% Constraints (including variable bounds):
%   g(1): (0.5*(x(1)+x(2))^2+150)-267.42                <= 0;
%   g(2): -(0.5*(x(1)+x(2))^2+150)                      <= 0;
%   h(1): 30*x(1)-6*x(1)^2-(0.5*(x(1)+x(2))^2+150)+250   = 0;
%   h(2): 20*x(2)-12*x(2)^2-(0.5*(x(1)+x(2))^2+150)+300  = 0;
%         0 <= x(1) <= 9.422;
%         0 <= x(2) <= 5.903;
%   
% Problem Properties:
%   n  = 2;
%   #g = 2;
%   #h = 2;  
% ------------------------------------------------------------------------------ 
y = 0.5*((x(1)+x(2))^2)+150;  
end