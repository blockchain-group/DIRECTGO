function y = P3b(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   P3b.m
%
% Original source: 
% - Birgin, E. G., Floudas, C. A., Martínez, J. M. :Global minimization usingan 
%   Augmented Lagrangian method with variable lower-level constraints. Math. 
%   Program. Ser. A 125(1), 139–162 (2010)
%
% Globally optimal solution:
%   f* = -0.38881
%   x* = (3.03556750332836, 5.09726343577258) 
%
% Constraints (including variable bounds):
%   g(1): x(1)^(1/2)+x(1)^(1/2)-4 <= 0;
%         10^(-5) <= x(1) <= 16;
%         10^(-5) <= x(1) <= 16;
%   
% Problem Properties:
%   n  = 2;
%   #g = 1;
%   #h = 0;  
% ------------------------------------------------------------------------------ 
k1 = 0.09755988;
k2 = 0.99*k1;
k3 = 0.03919080;
k4 = 0.9*k3;
y = -(((k1*x(1))/((1+k1*x(1))*(1+k3*x(1))*(1+k4*x(2))))+((k2*x(2))/((1+k1*x(1))*(1+k2*x(2))*(1+k4*x(2))))); 
end