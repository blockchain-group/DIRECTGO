function y = Horst1(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Horst1.m
%
% Original source: 
% - Horst, R., Pardalos, P.M., Thoai, N.V. (1995). Introduction to  
%   Global Optimization. Nonconvex Optimization and Its Application. 
%   Kluwer, Dordrecht  
%
% Globally optimal solution:
%   f* = -1.0625
%   x* = (0.75, 2) 
%
% Constraints (including variable bounds):
%   g(1): x(1)-4*x(2)    <= 1;
%   g(2): x(1)+x(2)      <= 4;
%   g(3): -4*x(1)+2*x(2) <= 1;
%         0 <= x(1) <= 3;
%         0 <= x(2) <= 2;
%   
% Problem Properties:
%   n  = 2;
%   #g = 3;
%   #h = 0;  
% ------------------------------------------------------------------------------ 
y = -x(1)^2-4*x(2)^2+4*x(1)*x(2)+2*x(1)+4*x(2);
end