function y = Horst5(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Horst5.m
%
% Original source: 
% - Horst, R., Pardalos, P.M., Thoai, N.V. (1995). Introduction to  
%   Global Optimization. Nonconvex Optimization and Its Application. 
%   Kluwer, Dordrecht  
%
% Globally optimal solution:
%   f* = -3.7220393738285287
%   x* = (1.2, 0, 0.8)
%
% Constraints (including variable bounds):
%   g(1): x(3)                 <= 3;
%   g(2): -2*x(1)-2*x(2)+x(3)  <= 1;
%   g(3): x(1)+x(2)-(1/4)*x(3) <= 1;
%   g(4): x(1)+x(2)+x(3)       <= 2;
%         0 <= x(1) <= 1.2;
%         0 <= x(2) <= 1.2;
%         0 <= x(3) <= 1.7;
%   
% Problem Properties:
%   n  = 3;
%   #g = 4;
%   #h = 0;  
% -------------------------------------------------------------------------
if nargin == 0
    y.nx = 3;
    y.ng = 4;
    y.nh = 0;
    y.xl = @(i) 0;
    xu = [1.2, 1.2, 1.7];
    y.xu = @(i) xu(i);
    y.fmin = @(i) -3.7220393738285287;
    xmin = [1.2, 0, 0.8];
    y.xmin = @(i) xmin(i);
    y.confun = @(i) Horst5c(i);
    return
end
y = (-(abs(x(1) + (1/2)*x(2) + (2/3)*x(3)))^(3/2)) - x(1)^2;
end

function [c, ceq] = Horst5c( x )
c(1) = x(1) + x(2) + x(3) - 2;
c(2) = x(1) + x(2) - (1/4)*x(3) - 1;
c(3) = -2*x(1) - 2*x(2) + x(3) - 1;
c(4) = x(3) - 3;
ceq = [];
end