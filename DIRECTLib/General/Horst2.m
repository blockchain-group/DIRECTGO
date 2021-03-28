function y = Horst2(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Horst2.m
%
% Original source: 
% - Horst, R., Pardalos, P.M., Thoai, N.V. (1995). Introduction to  
%   Global Optimization. Nonconvex Optimization and Its Application. 
%   Kluwer, Dordrecht  
%
% Globally optimal solution:
%   f* = -6.8995
%   x* = (2.5, 0.75) 
%
% Constraints (including variable bounds):
%   g(1): -x(1)+x(2)  <= 1;
%   g(2): x(1)-2*x(2) <= 1;
%   g(3): x(1)+2*x(2) <= 4;
%         0 <= x(1) <= 2.5;
%         0 <= x(2) <= 2;
%   
% Problem Properties:
%   n  = 2;
%   #g = 3;
%   #h = 0;  
% ------------------------------------------------------------------------------ 
y = -x(1)^2-x(2)^(3/2);
end