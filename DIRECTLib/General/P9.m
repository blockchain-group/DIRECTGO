function y = P9(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   P9.m
%
% Original source: 
% - Birgin, E. G., Floudas, C. A., Martínez, J. M. :Global minimization usingan 
%   Augmented Lagrangian method with variable lower-level constraints. Math. 
%   Program. Ser. A 125(1), 139–162 (2010)
%
% Globally optimal solution:
%   f* = -13.401536871200
%   x* = (0.166666666666, 2, 4) 
%
% Constraints (including variable bounds):
%   g(1): -4*x(1)+(4/3)*x(2)-6        <= 0;
%   g(2): -x(2)+(1/2)*x(3)-2          <= 0;
%   g(3): -x(1)+(1/3)*x(2)-2          <= 0;
%   g(4): x(1)+2*(-x(1)+(1/3)*x(2))-4 <= 0;
%   g(5): x(2)+(-x(2)+(1/2)*x(3))- 4  <= 0;
%   g(6): x(3)+(-4*x(1)+(4/3)*x(2))-6 <= 0;
%   g(7): -(-4*x(1)+(4/3)*x(2))       <= 0;
%   g(8): -(-x(2)+(1/2)*x(3))         <= 0;
%   g(9): -(-x(1)+(1/3)*x(2))         <= 0;
%         10^(-5) <= x(1) <= 3;
%         10^(-5) <= x(2) <= 4;
%         10^(-5) <= x(3) <= 4;
%   
% Problem Properties:
%   n  = 3;
%   #g = 9;
%   #h = 0;  
% ------------------------------------------------------------------------------ 
y = x(1)^(0.6)+x(2)^(0.6)+x(3)^(0.4)-(3/2)*x(3)+2*x(1)-(17/3)*x(2); 
end