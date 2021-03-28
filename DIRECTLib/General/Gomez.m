function y = Gomez(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Gomez.m
%
% Original source: 
% - Birgin, E. G., Floudas, C. A., Martínez, J. M. :Global minimization usingan 
%   Augmented Lagrangian method with variable lower-level constraints. Math. 
%   Program. Ser. A 125(1), 139–162 (2010)
%
% Globally optimal solution:
%   f* = -0.971475938409664
%   x* = (0.10928197498268, -0.62375796787744) 
%
% Constraints (including variable bounds):
%   g(1): -sin(4*3.14*x(1)) + 2*sin(2*3.14*x(2))^2 <= 0;
%         -1 <= x(1) <= 1;
%         -1 <= x(2) <= 1;
%   
% Problem Properties:
%   n  = 2;
%   #g = 1;
%   #h = 0;  
% ------------------------------------------------------------------------------ 
y = (4-2.1*x(1)^2+(x(1)^4)/3)*x(1)^2+x(1)*x(2)+(-4+4*x(2)^2)*x(2)^2;
end