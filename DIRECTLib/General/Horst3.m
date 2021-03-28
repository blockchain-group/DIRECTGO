function y = Horst3(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Horst3.m
%
% Original source: 
% - Horst, R., Pardalos, P.M., Thoai, N.V. (1995). Introduction to  
%   Global Optimization. Nonconvex Optimization and Its Application. 
%   Kluwer, Dordrecht  
%
% Globally optimal solution:
%   f* = -(4/9)
%   x* = (0, 0) 
%
% Constraints (including variable bounds):
%   g(1): x(1)+(1/10)*x(2) <= 1;
%   g(2): x(1)+x(2)        <= (3/2);
%   g(3): -2*x(1)+x(2)     <= 1;
%         0 <= x(1) <= 1;
%         0 <= x(2) <= 1.5;
%   
% Problem Properties:
%   n  = 2;
%   #g = 3;
%   #h = 0;  
% ------------------------------------------------------------------------------ 
y = -x(1)^2+(4/3)*x(1)+((log(1+x(2)))/(log(exp(1))))-(4/9);
end