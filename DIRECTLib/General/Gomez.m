function y = Gomez(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Gomez.m
%
% Original source: 
% - Birgin, E. G., Floudas, C. A., Martínez, J. M. :Global minimization 
%   usingan Augmented Lagrangian method with variable lower-level 
%   constraints. Math. Program. Ser. A 125(1), 139–162 (2010)
%
% Globally optimal solution:
%   f* = -0.9714759185876088
%   x* = (0.1092819737821463, -0.6237579513444821) 
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
% -------------------------------------------------------------------------
if nargin == 0
    y.nx = 2;
    y.ng = 1;
    y.nh = 0;
    y.xl = @(i) -1;
    y.xu = @(i) 1;
    y.fmin = @(i) -0.9714759185876088;
    xmin = [0.1092819737821463, -0.6237579513444821];
    y.xmin = @(i) xmin(i);
    y.confun = @(i) Gomezc(i);
    return
end
y = (4 - 2.1*x(1)^2 + (x(1)^4)/3)*x(1)^2 + x(1)*x(2) + (-4 +...
    4*x(2)^2)*x(2)^2;
end

function [c, ceq] = Gomezc(x)
c   = -sin(4*3.14*x(1)) + 2*sin(2*3.14*x(2))^2; 
ceq = [];
end

