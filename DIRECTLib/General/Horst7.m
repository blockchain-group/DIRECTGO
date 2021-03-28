function y = Horst7(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Horst7.m
%
% Original source: 
% - Horst, R., Pardalos, P.M., Thoai, N.V. (1995). Introduction to  
%   Global Optimization. Nonconvex Optimization and Its Application. 
%   Kluwer, Dordrecht  
%
% Globally optimal solution:
%   f* = -52.8774169979695
%   x* = (6, 0, 3)
%
% Constraints (including variable bounds):
%   g(1): x(3)                  <= 3;
%   g(2): -2*x(1)-4*x(2)-2*x(3) <= -1;
%   g(3): x(1)+2*x(2)           <= 6;
%   g(4): -x(1)-x(2)+(1/2)*x(3) <= 1;
%         0 <= x(1) <= 6;
%         0 <= x(2) <= 3;
%         0 <= x(3) <= 3;
%   
% Problem Properties:
%   n  = 3;
%   #g = 4;
%   #h = 0;  
% ------------------------------------------------------------------------------ 
y = -((x(1)+(1/2)*x(3)-2)^2)-(abs(x(1)+(1/2)*x(2)+(2/3)*x(3)))^(3/2);
end