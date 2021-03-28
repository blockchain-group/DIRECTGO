function y = P3a(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   P3a.m
%
% Original source: 
% - Birgin, E. G., Floudas, C. A., Martínez, J. M. :Global minimization usingan 
%   Augmented Lagrangian method with variable lower-level constraints. Math. 
%   Program. Ser. A 125(1), 139–162 (2010)
%
% Globally optimal solution:
%   f* = -0.38881
%   x* = (0.772, 0.517, 0.204, 0.388, 3.036, 5.097) 
%
% Constraints (including variable bounds):
%   g(1): x(5)^(1/2)+x(6)^(1/2)-4                 <= 0;
%   h(1): x(1)+0.09755988*x(1)*x(5)-1              = 0;
%   h(2): x(2)-x(1)+0.0965842812*x(2)*x(6)         = 0;
%   h(3): x(3)+x(1)+0.03919080*x(3)*x(5)-1         = 0;
%   h(4): x(4)-x(3)+x(2)-x(1)+0.03527172*x(4)*x(6) = 0;
%         0       <= x(1) <= 1;
%         0       <= x(2) <= 1;
%         0       <= x(3) <= 1;
%         0       <= x(4) <= 1;
%         10^(-5) <= x(5) <= 16;
%         10^(-5) <= x(6) <= 16;
%   
% Problem Properties:
%   n  = 6;
%   #g = 1;
%   #h = 4;  
% ------------------------------------------------------------------------------ 
y = -x(4);
end