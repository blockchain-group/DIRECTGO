function y = P7(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   P7.m
%
% Original source: 
% - Birgin, E. G., Floudas, C. A., Martínez, J. M. :Global minimization usingan 
%   Augmented Lagrangian method with variable lower-level constraints. Math. 
%   Program. Ser. A 125(1), 139–162 (2010)
%
% Globally optimal solution:
%   f* = -2.8284
%   x* = (-1.41421356187310, -1.41421356187309) 
%
% Constraints (including variable bounds):
%   g(1): -x(1)+x(2)-1     <= 0;
%   g(2): x(1)-x(2)-1      <= 0;
%   g(3): -x(1)^2-x(2)^2+1 <= 0;
%   g(4): x(1)^2+x(2)^2-4  <= 0;
%         -2 <= x(1) <= 2;
%         -2 <= x(2) <= 2;
%   
% Problem Properties:
%   n  = 2;
%   #g = 4;
%   #h = 0;  
% ------------------------------------------------------------------------------ 
y = x(1)+x(2); 
end