function y = P1(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   P1.m
%
% Original source: 
% - Birgin, E. G., Floudas, C. A., Martínez, J. M. :Global minimization usingan 
%   Augmented Lagrangian method with variable lower-level constraints. Math. 
%   Program. Ser. A 125(1), 139–162 (2010)
%
% Globally optimal solution:
%   f* = 0.0292840682010227
%   x* = (1.11659549992298, 1.22038649820888, 1.53782370706329,
%         1.97284238054318, 1.79106937126107) 
%
% Constraints (including variable bounds):
%   g(1): x(1)+x(2)^2+x(3)^3-3*sqrt(2)-2 = 0;
%   g(2): x(2)-x(3)^2+x(4)-2*sqrt(2)+2   = 0;
%   g(3): x(1)*x(5)-2                    = 0;
%         -5 <= x(1) <= 5;
%         -5 <= x(2) <= 5;
%         -5 <= x(3) <= 5;
%         -5 <= x(4) <= 5;
%         -5 <= x(5) <= 5;
%   
% Problem Properties:
%   n  = 5;
%   #g = 0;
%   #h = 3;  
% ------------------------------------------------------------------------------ 
y = (x(1)-1)^2+(x(1)-x(2))^2+(x(2)-x(3))^3+(x(3)-x(4))^4+(x(4)-x(5))^4;
end