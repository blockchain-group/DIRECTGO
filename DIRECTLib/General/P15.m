function y = P15(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   P15.m
%
% Original source: 
% - Birgin, E. G., Floudas, C. A., Martínez, J. M. :Global minimization usingan 
%   Augmented Lagrangian method with variable lower-level constraints. Math. 
%   Program. Ser. A 125(1), 139–162 (2010)
%
% Globally optimal solution:
%   f* = 0
%   x* = (10.416668333333, 31.250001666666, 8.3333333333333) 
%
% Constraints (including variable bounds):
%   h(1): x(1)+x(2)+x(3)-50             = 0;
%   h(2): x(2)/x(1)-3                   = 0;
%   h(3): x(3)^2/(x(1)*x(2)^3)-0.000169 = 0;
%         10^(-5) <= x(1) <= 12.5;
%         10^(-5) <= x(2) <= 37.5;
%         0       <= x(3) <= 50;
%   
% Problem Properties:
%   n  = 3;
%   #g = 0;
%   #h = 3;  
% ------------------------------------------------------------------------------ 
y = 0; 
end