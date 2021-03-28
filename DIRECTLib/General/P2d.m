function y = P2d(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   P2d_mod.m
%
% Original source: 
% - Birgin, E. G., Floudas, C. A., Martínez, J. M. :Global minimization usingan 
%   Augmented Lagrangian method with variable lower-level constraints. Math. 
%   Program. Ser. A 125(1), 139–162 (2010)
%
% Globally optimal solution:
%   f* = -400
%   x* = (0, 100, 0, 100, 1) 
%
% Constraints (including variable bounds):
%   g(1): x(5)*x(1)+2*x(3)-2.5*(x(1) + x(3))                   <= 0;
%   g(2): x(5)*x(2)+2*x(4)-1.5*(x(2)+x(4))                     <= 0;
%   g(3): x(3)+x(4)-300                                        <= 0;
%   g(4): -(x(3)+x(4))                                         <= 0;
%   g(5): x(2)+x(4) - 200                                      <= 0;
%   g(6): -(x(2)+x(4))                                         <= 0;
%   g(7): x(1) + x(3) - 100                                    <= 0;
%   g(8): -(x(1) + x(3))                                       <= 0;
%   g(9): ((x(1)*x(5)+x(2)*x(5)-x(1)-x(2))/2) -300             <= 0;
%   g(10): -((x(1)*x(5)+x(2)*x(5)-x(1)-x(2))/2)                <= 0;
%   g(11): (x(1)+x(2)-((x(1)*x(5)+x(2)*x(5)-x(1)-x(2))/2))-300 <= 0;
%   g(12): -((x(1)+x(2)-((x(1)*x(5)+x(2)*x(5)-x(1)-x(2))/2)))  <= 0;
%         0 <= x(1) <= 100;
%         0 <= x(2) <= 200;
%         0 <= x(3) <= 100;
%         0 <= x(4) <= 200;
%         1 <= x(5) <= 3;
%   
% Problem Properties:
%   n  = 5;
%   #g = 12;
%   #h = 0;  
% ------------------------------------------------------------------------------ 
y = -9*(x(1)+x(3))-15*(x(2)+x(4))+6*((x(5)*(x(1)+x(2))-x(1)-x(2))/2)+16*((x(1)+x(2)-((x(5)*(x(1)+x(2))-x(1)-x(2))/2)))+10*(x(3)+x(4)); 
end