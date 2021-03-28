function y = P2b(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   P2b.m
%
% Original source: 
% - Birgin, E. G., Floudas, C. A., Martínez, J. M. :Global minimization usingan 
%   Augmented Lagrangian method with variable lower-level constraints. Math. 
%   Program. Ser. A 125(1), 139–162 (2010)
%
% Globally optimal solution:
%   f* = -400
%   x* = (0, 100, 1, 0, 100) 
%
% Constraints (including variable bounds):
%   g(1): x(4)+x(1)-600                                        <= 0;
%   g(2): -(x(4)+x(1))                                         <= 0;
%   g(3): x(5)+x(2)-200                                        <= 0;
%   g(4): -(x(5)+x(2))                                         <= 0;
%   g(5): x(3)*x(5) + 2*x(2) - 1.5*(x(5)+x(2))                 <= 0;
%   g(6): x(3)*x(4) + 2*x(1) - 2.5*(x(4)+x(1))                 <= 0;
%   g(7): ((x(3)*x(4)+x(3)*x(5)-x(4)-x(5))/2)-500              <= 0;
%   g(8): -((x(3)*x(4)+x(3)*x(5)-x(4)-x(5))/2)                 <= 0;
%   g(9): x(4)+x(5) - ((x(3)*x(4)+x(3)*x(5)-x(4)-x(5))/2) -500 <= 0;
%   g(10): -(x(4)+x(5) - ((x(3)*x(4)+x(3)*x(5)-x(4)-x(5))/2))  <= 0;
%         0 <= x(1) <= 500;
%         0 <= x(2) <= 500;
%         0 <= x(3) <= 500;
%         0 <= x(4) <= 500;
%         0 <= x(5) <= 500;
%   
% Problem Properties:
%   n  = 5;
%   #g = 10;
%   #h = 0;  
% ------------------------------------------------------------------------------ 
y = -9*(x(4)+x(1))-15*(x(5)+x(2))+6*(((x(3)*x(4)+x(3)*x(5)-x(4)-x(5))/2))+16*(x(4)+x(5)-((x(3)*x(4)+x(3)*x(5)-x(4)-x(5))/2))+10*(x(1)+x(2)); 
end