function y = P13(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   P13.m
%
% Original source: 
% - Birgin, E. G., Floudas, C. A., Martínez, J. M. :Global minimization usingan 
%   Augmented Lagrangian method with variable lower-level constraints. Math. 
%   Program. Ser. A 125(1), 139–162 (2010)
%
% Globally optimal solution:
%   f* = 189.35
%   x* = (0.717536188588019) 
%
% Constraints (including variable bounds):
%   h(1): 600*x(1)-50*x(3)-x(1)*x(3)+5000 = 0;
%   h(2): 600*x(2)+50*x(3)-15000          = 0;
%         10^(-5) <= x(1) <= 34;
%         10^(-5) <= x(2) <= 17;
%         100     <= x(3) <= 300;
%   
% Problem Properties:
%   n  = 3;
%   #g = 0;
%   #h = 2;  
% ------------------------------------------------------------------------------ 
y = 35*x(1)^(0.6)+35*x(2)^(0.6); 
end